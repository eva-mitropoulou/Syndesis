from __future__ import annotations

import hashlib
import importlib.metadata
import json
import os
import warnings as py_warnings
from dataclasses import dataclass
from pathlib import Path
import shutil
import subprocess
from typing import Any

import MDAnalysis as mda
import numpy as np
import pandas as pd
import prolif as plf
from openmm.app import PDBFile
from pdbfixer import PDBFixer

from syndesis.stage5.residue_mapping import residue_role

_PROTEIN_MOL_CACHE: dict[str, plf.Molecule] = {}


@dataclass(frozen=True)
class Atom:
    index: int
    name: str
    residue_name: str
    chain_id: str
    residue_number: int | None
    element: str
    xyz: np.ndarray


def interaction_config_hash(config: dict[str, Any]) -> str:
    payload = json.dumps(config.get("interactions", {}), sort_keys=True)
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()[:16]


def engine_metadata() -> tuple[str, str]:
    try:
        return "prolif", importlib.metadata.version("prolif")
    except importlib.metadata.PackageNotFoundError:
        raise RuntimeError("Stage 5 requires ProLIF. Install `prolif` and rerun; no fallback engine is allowed.")


def _element_from_atom_line(line: str) -> str:
    element = line[76:78].strip()
    if element:
        return element.upper()
    name = line[12:16].strip()
    letters = "".join(ch for ch in name if ch.isalpha())
    return (letters[:1] or "C").upper()


def read_atoms(path: str | Path) -> list[Atom]:
    atoms: list[Atom] = []
    for raw in Path(path).read_text(encoding="utf-8", errors="ignore").splitlines():
        if not raw.startswith(("ATOM", "HETATM")):
            continue
        try:
            xyz = np.array([float(raw[30:38]), float(raw[38:46]), float(raw[46:54])], dtype=float)
        except ValueError:
            continue
        resnum_text = raw[22:26].strip()
        try:
            resnum = int(resnum_text)
        except ValueError:
            resnum = None
        atoms.append(
            Atom(
                index=len(atoms),
                name=raw[12:16].strip(),
                residue_name=raw[17:20].strip(),
                chain_id=raw[21:22].strip(),
                residue_number=resnum,
                element=_element_from_atom_line(raw),
                xyz=xyz,
            )
        )
    return atoms


def has_hydrogens(path: str | Path) -> bool:
    return any(atom.element == "H" for atom in read_atoms(path))


def _openbabel_environment(executable: str) -> dict[str, str] | None:
    plugin_dir = Path(executable).parent.parent / "lib" / "openbabel" / "3.1.0"
    return {**os.environ, "BABEL_LIBDIR": str(plugin_dir)} if plugin_dir.exists() else None


def _autodock_element(line: str) -> str:
    atom_type = line[77:].strip().split()[0] if len(line) > 77 and line[77:].strip() else ""
    mapping = {
        "A": "C",
        "C": "C",
        "N": "N",
        "NA": "N",
        "OA": "O",
        "O": "O",
        "SA": "S",
        "S": "S",
        "HD": "H",
        "H": "H",
        "F": "F",
        "CL": "CL",
        "Cl": "CL",
        "BR": "BR",
        "Br": "BR",
        "I": "I",
    }
    if atom_type in mapping:
        return mapping[atom_type]
    atom_name = line[12:16].strip()
    letters = "".join(ch for ch in atom_name if ch.isalpha())
    if len(letters) >= 2 and letters[:2].title() in {"Cl", "Br"}:
        return letters[:2].upper()
    return (letters[:1] or "C").upper()


def prepare_ligand_for_prolif(ligand_file: str | Path, out_dir: str | Path) -> Path:
    source = Path(ligand_file)
    target_dir = Path(out_dir)
    target_dir.mkdir(parents=True, exist_ok=True)
    source_key = hashlib.sha256(str(source.resolve()).encode("utf-8")).hexdigest()[:10]
    target = target_dir / f"{source.parent.name}__{source.stem}__{source_key}.prolif_h.pdb"
    if target.exists() and target.stat().st_mtime >= source.stat().st_mtime:
        return target
    obabel = shutil.which("obabel") or str(Path(".tools/stage11_amber/bin/obabel").resolve())
    if not Path(obabel).exists():
        raise RuntimeError("OpenBabel `obabel` is required to prepare explicit-hydrogen ligands for ProLIF.")
    input_format = "pdbqt" if source.suffix.lower() == ".pdbqt" else "pdb" if source.suffix.lower() == ".pdb" else None
    if input_format is None:
        input_format = source.suffix.lower().lstrip(".")
    command = [
        obabel,
        f"-i{input_format}",
        str(source),
        "-opdb",
        "-O",
        str(target),
        "-h",
        "-p",
        "7.4",
    ]
    proc = subprocess.run(command, text=True, capture_output=True, check=False, env=_openbabel_environment(obabel))
    if proc.returncode != 0:
        raise RuntimeError(f"OpenBabel failed to add ligand hydrogens for ProLIF: {source}\n{proc.stderr}")
    if not target.exists() or target.stat().st_size == 0:
        raise RuntimeError(f"OpenBabel produced no ProLIF ligand file: {target}")
    if not has_hydrogens(target):
        raise RuntimeError(f"Prepared ProLIF ligand lacks explicit hydrogens: {target}")
    return target


def prepare_protein_for_prolif(protein_file: str | Path, out_dir: str | Path, ph: float = 7.4) -> Path:
    source = Path(protein_file)
    target_dir = Path(out_dir)
    target_dir.mkdir(parents=True, exist_ok=True)
    target = target_dir / f"{source.parent.name}__{source.stem}.h.pdb"
    if target.exists() and target.stat().st_mtime >= source.stat().st_mtime:
        return target
    # Docked-pose interaction analysis must use precisely the receptor heavy atoms
    # supplied to docking. When the Stage 3 PDBQT sibling exists, it is the
    # authoritative Open Babel-protonated representation. Converting it back to
    # PDB retains its heavy-atom set for ProLIF; PDBFixer is deliberately not
    # allowed to repair missing heavy atoms in this branch.
    docking_pdbqt = source if source.suffix.lower() == ".pdbqt" else source.with_suffix(".pdbqt")
    if docking_pdbqt.exists():
        obabel = shutil.which("obabel") or str(Path(".tools/stage11_amber/bin/obabel").resolve())
        if not Path(obabel).exists():
            raise RuntimeError("OpenBabel `obabel` is required to derive the ProLIF receptor from docking PDBQT.")
        command = [obabel, "-ipdbqt", str(docking_pdbqt), "-opdb", "-O", str(target), "-h"]
        proc = subprocess.run(command, text=True, capture_output=True, check=False, env=_openbabel_environment(obabel))
        if proc.returncode != 0 or not target.exists() or target.stat().st_size == 0:
            raise RuntimeError(f"OpenBabel failed to derive ProLIF receptor from docking PDBQT: {docking_pdbqt}\n{proc.stderr}")
        if not has_hydrogens(target):
            raise RuntimeError(f"Docking-derived ProLIF receptor lacks explicit hydrogens: {target}")
        return target
    fixer = PDBFixer(filename=str(source))
    try:
        fixer.addMissingHydrogens(ph)
    except ValueError as exc:
        raise RuntimeError(
            "PDBFixer could not add hydrogens without repairing missing protein heavy atoms; "
            f"provide the docking PDBQT sibling or a complete receptor: {source}"
        ) from exc
    with target.open("w", encoding="utf-8") as handle:
        PDBFile.writeFile(fixer.topology, fixer.positions, handle, keepIds=True)
    if not has_hydrogens(target):
        raise RuntimeError(f"Prepared ProLIF protein lacks explicit hydrogens: {target}")
    return target


def _parse_residue_id(value: Any) -> tuple[str, int | None, str]:
    text = str(value)
    residue_part, _, chain = text.partition(".")
    letters = "".join(ch for ch in residue_part if ch.isalpha())
    digits = "".join(ch for ch in residue_part if ch.isdigit() or ch == "-")
    try:
        resnum = int(digits)
    except ValueError:
        resnum = None
    return letters.upper(), resnum, chain


def _metadata_values(payload: Any) -> list[dict[str, Any]]:
    if isinstance(payload, tuple):
        return [item for item in payload if isinstance(item, dict)]
    if isinstance(payload, list):
        return [item for item in payload if isinstance(item, dict)]
    if isinstance(payload, dict):
        return [payload]
    return []


def _prolif_parameters(config: dict[str, Any]) -> dict[str, dict[str, Any]]:
    cutoffs = config.get("interactions", {}).get("distance_cutoffs", {})
    params: dict[str, dict[str, Any]] = {}
    if "hydrophobic" in cutoffs:
        params["Hydrophobic"] = {"distance": float(cutoffs["hydrophobic"])}
    if "hydrogen_bond" in cutoffs:
        params["ImplicitHBDonor"] = {
            "distance": float(cutoffs["hydrogen_bond"]),
            "ignore_geometry_checks": True,
        }
        params["ImplicitHBAcceptor"] = {
            "distance": float(cutoffs["hydrogen_bond"]),
            "ignore_geometry_checks": True,
        }
    if "ionic" in cutoffs:
        params["Anionic"] = {"distance": float(cutoffs["ionic"])}
        params["Cationic"] = {"distance": float(cutoffs["ionic"])}
        params["CationPi"] = {"distance": float(cutoffs["ionic"])}
        params["PiCation"] = {"distance": float(cutoffs["ionic"])}
    if "vdw" in cutoffs:
        params["VdWContact"] = {"tolerance": float(cutoffs["vdw"]) - 3.0}
    return params


def compute_interactions(
    protein_file: str | Path,
    ligand_file: str | Path,
    residue_map: pd.DataFrame,
    config: dict[str, Any],
) -> tuple[pd.DataFrame, dict[str, Any]]:
    engine, version = engine_metadata()
    warnings: list[str] = []
    if not has_hydrogens(protein_file):
        warnings.append("protein_hydrogens_missing_or_not_explicit")
    if not has_hydrogens(ligand_file):
        warnings.append("ligand_hydrogens_missing_or_not_explicit")
    if config.get("interactions", {}).get("require_hydrogens", False) and warnings:
        raise RuntimeError(f"Explicit hydrogens are required for ProLIF interactions: {warnings}")
    protein_key = str(Path(protein_file).resolve())
    ligand_universe = mda.Universe(str(ligand_file))
    with py_warnings.catch_warnings():
        py_warnings.filterwarnings("ignore", message="No hydrogen atom found.*")
        py_warnings.filterwarnings("ignore", message="No `bonds` attribute.*")
        if protein_key not in _PROTEIN_MOL_CACHE:
            protein_universe = mda.Universe(str(protein_file))
            _PROTEIN_MOL_CACHE[protein_key] = plf.Molecule.from_mda(protein_universe, inferrer=None, force=True)
        protein = _PROTEIN_MOL_CACHE[protein_key]
        ligand = plf.Molecule.from_mda(ligand_universe, force=True)
    requested = config.get("interactions", {}).get("enabled_interactions") or None
    available = set(plf.Fingerprint.list_available())
    interactions = [name for name in requested if name in available] if requested else None
    unsupported = sorted(set(requested or []) - available)
    if unsupported:
        raise RuntimeError(f"Unsupported ProLIF interactions requested: {unsupported}")
    fp = plf.Fingerprint(
        interactions=interactions,
        parameters=_prolif_parameters(config),
        implicit_hydrogens=True,
    )
    ifp = fp.generate(ligand, protein, metadata=True)
    residue_lookup = {
        int(float(row["auth_seq_id"])): row
        for _, row in residue_map.iterrows()
        if pd.notna(row.get("auth_seq_id"))
    }
    rows: list[dict[str, Any]] = []
    for (_lig_residue, protein_residue), interaction_payload in ifp.items():
        resname, resnum, _chain = _parse_residue_id(protein_residue)
        mapped = residue_lookup.get(resnum) if resnum is not None else None
        uniprot = mapped.get("uniprot_residue_number") if mapped is not None else resnum
        klifs = mapped.get("klifs_position") if mapped is not None else None
        role = mapped.get("residue_role") if mapped is not None else residue_role(uniprot)
        for interaction_type, metadata_payload in interaction_payload.items():
            for metadata in _metadata_values(metadata_payload):
                indices = metadata.get("indices", {})
                parents = metadata.get("parent_indices", {})
                distance = metadata.get("distance")
                rows.append(
                    {
                        "residue_name": resname,
                        "auth_seq_id": resnum,
                        "uniprot_residue_number": uniprot,
                        "klifs_position": klifs,
                        "residue_role": role,
                        "interaction_type": interaction_type,
                        "present": True,
                        "atom_indices_ligand_json": json.dumps(list(indices.get("ligand", parents.get("ligand", [])))),
                        "atom_indices_protein_json": json.dumps(list(indices.get("protein", parents.get("protein", [])))),
                        "distance": None if distance is None else round(float(distance), 3),
                        "angle": metadata.get("angle"),
                        "interaction_confidence": "medium" if warnings else "high",
                        "warnings_json": json.dumps(warnings),
                    }
                )
    meta = {
        "interaction_engine": engine,
        "interaction_engine_version": version,
        "interaction_config_hash": interaction_config_hash(config),
        "warnings": warnings,
    }
    return pd.DataFrame(rows), meta


def fingerprint_from_interactions(interactions: pd.DataFrame) -> tuple[str, str, set[str]]:
    if interactions.empty:
        return "", "[]", set()
    bits = []
    for row in interactions.itertuples(index=False):
        residue = row.uniprot_residue_number if pd.notna(row.uniprot_residue_number) else row.auth_seq_id
        try:
            residue_text = str(int(float(residue)))
        except (TypeError, ValueError):
            residue_text = "NA"
        bits.append(f"{row.residue_role}:{residue_text}:{row.interaction_type}")
    bits = sorted(bits)
    return "|".join(bits), json.dumps(bits), set(bits)


def tanimoto(a: set[str], b: set[str]) -> float | None:
    if not a and not b:
        return None
    union = a | b
    if not union:
        return None
    return float(len(a & b) / len(union))
