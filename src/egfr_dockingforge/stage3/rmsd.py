from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from egfr_dockingforge.common.io import write_table
from egfr_dockingforge.stage3.schemas import RMSD_COLUMNS


def heavy_atom_coords(path: str | Path) -> np.ndarray:
    coords = []
    with Path(path).open("r", encoding="utf-8") as handle:
        for line in handle:
            if not line.startswith(("ATOM", "HETATM")):
                continue
            element = line[76:78].strip() or line[12:16].strip()[0]
            atom_type = line[77:].strip().split()[0] if len(line) > 77 and line[77:].strip() else element
            if element.upper() in {"H", "HD"} or atom_type.upper() in {"H", "HD"}:
                continue
            coords.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
    return np.asarray(coords, dtype=float)


def mapped_rmsd(a: np.ndarray, b: np.ndarray) -> float:
    """Naive positional RMSD assuming identical atom ordering (no symmetry)."""
    if len(a) != len(b):
        raise ValueError(f"Atom count mismatch: {len(a)} vs {len(b)}")
    if len(a) == 0:
        raise ValueError("No heavy atoms available for RMSD.")
    return float(np.sqrt(((a - b) ** 2).sum(axis=1).mean()))


def _load_mol_with_coords(path: str | Path, smiles: str | None):
    """Load a docked/reference ligand as an RDKit mol with 3D coordinates.

    Docked poses are PDB/PDBQT with unreliable bond orders, so when a reference
    SMILES is available we assign bond orders from it (this also fixes element
    perception and lets us use the canonical atom-matching machinery). Returns
    ``None`` if a chemically valid molecule cannot be built.
    """
    from rdkit import Chem
    from rdkit.Chem import AllChem

    path = str(path)
    mol = None
    if path.lower().endswith(".pdbqt"):
        # RDKit cannot read PDBQT directly; strip the trailing PDBQT-only
        # columns so the ATOM/HETATM records parse as PDB.
        pdb_block_lines = []
        for line in Path(path).read_text(encoding="utf-8", errors="ignore").splitlines():
            if line.startswith(("ATOM", "HETATM")):
                pdb_block_lines.append(line[:66])
            elif line.startswith(("END", "TER")):
                pdb_block_lines.append(line[:6])
        mol = Chem.MolFromPDBBlock("\n".join(pdb_block_lines), removeHs=True, sanitize=False)
    else:
        mol = Chem.MolFromPDBFile(path, removeHs=True, sanitize=False)
    if mol is None:
        return None
    if smiles:
        template = Chem.MolFromSmiles(smiles)
        if template is not None:
            template = Chem.RemoveHs(template)
            try:
                mol = AllChem.AssignBondOrdersFromTemplate(template, mol)
            except Exception:
                # Fall back to the raw perceived connectivity.
                pass
    try:
        Chem.SanitizeMol(mol)
    except Exception:
        return None
    return mol


def symmetry_corrected_rmsd(pose_path: str | Path, ref_path: str | Path, smiles: str | None) -> tuple[float, str, str]:
    """Symmetry-aware heavy-atom RMSD in the current (already-aligned) frame.

    Poses and reference are pre-placed in a common frame (the reference is
    transformed into the receptor frame upstream), so we must NOT re-superpose;
    we only need symmetry-aware atom matching. Returns (rmsd, method, status).
    """
    from rdkit.Chem import rdMolAlign

    pose = _load_mol_with_coords(pose_path, smiles)
    ref = _load_mol_with_coords(ref_path, smiles)
    if pose is None or ref is None:
        # Fall back to the naive positional RMSD on raw heavy-atom coordinates.
        value = mapped_rmsd(heavy_atom_coords(pose_path), heavy_atom_coords(ref_path))
        return value, "naive_atom_order_fallback", "fallback_no_molgraph"
    if pose.GetNumAtoms() != ref.GetNumAtoms():
        value = mapped_rmsd(heavy_atom_coords(pose_path), heavy_atom_coords(ref_path))
        return value, "naive_atom_order_fallback", "fallback_atom_count_mismatch"
    # CalcRMS enumerates symmetry-equivalent atom mappings and returns the best,
    # WITHOUT modifying coordinates (no superposition) -> in-place symmetry RMSD.
    try:
        value = float(rdMolAlign.CalcRMS(pose, ref))
        return value, "rdkit_symmetry_calcrms_inplace", "symmetry_matched"
    except Exception as exc:  # noqa: BLE001
        value = mapped_rmsd(heavy_atom_coords(pose_path), heavy_atom_coords(ref_path))
        return value, "naive_atom_order_fallback", f"fallback_calcrms_error:{str(exc)[:80]}"


def compute_pose_rmsd(poses: pd.DataFrame, tasks: pd.DataFrame, transforms: pd.DataFrame, config: dict[str, Any], paths: dict[str, Path], ligand_prep: pd.DataFrame | None = None) -> pd.DataFrame:
    rows = []
    task_lookup = tasks.set_index("docking_task_id")
    transform_lookup = transforms.set_index(["ligand_id", "target_receptor_id"])
    smiles_lookup: dict[str, str] = {}
    if ligand_prep is not None and not ligand_prep.empty:
        for rec in ligand_prep.to_dict("records"):
            smi = rec.get("prepared_smiles") or rec.get("standard_smiles")
            if isinstance(smi, str) and smi:
                smiles_lookup[str(rec.get("ligand_id"))] = smi
    for _, pose in poses.iterrows():
        task = task_lookup.loc[pose["docking_task_id"]]
        transform = transform_lookup.loc[(pose["ligand_id"], pose["target_receptor_id"])]
        status = "complete"
        warning = ""
        naive_value = None
        sym_value = None
        method = "symmetry_attempted_mapped_heavy_atom_fallback"
        mapping_status = "failed"
        smiles = smiles_lookup.get(str(pose["ligand_id"]))
        try:
            naive_value = mapped_rmsd(
                heavy_atom_coords(pose["pose_file"]),
                heavy_atom_coords(transform["transformed_reference_pose_file"]),
            )
        except Exception as exc:  # noqa: BLE001
            warning = str(exc)
        try:
            sym_value, method, mapping_status = symmetry_corrected_rmsd(
                pose["pose_file"], transform["transformed_reference_pose_file"], smiles
            )
        except Exception as exc:  # noqa: BLE001
            sym_value = naive_value
            method = "naive_atom_order_fallback"
            mapping_status = f"symmetry_failed:{str(exc)[:80]}"
            if not warning:
                warning = str(exc)
        if naive_value is None and sym_value is None:
            status = "failed"
        rows.append({
            "pose_id": pose["pose_id"],
            "docking_task_id": pose["docking_task_id"],
            "ligand_id": pose["ligand_id"],
            "target_receptor_id": pose["target_receptor_id"],
            "task_type": task["task_type"],
            "docking_engine": pose["docking_engine"],
            "pose_rank": pose["pose_rank"],
            "docking_score": pose["docking_score"],
            "rmsd_heavy_atom": naive_value,
            "rmsd_symmetry_corrected": sym_value if sym_value is not None else naive_value,
            "rmsd_method": method,
            "atom_mapping_status": mapping_status,
            "atom_mapping_warning": warning,
            "reference_pose_file": transform["transformed_reference_pose_file"],
            "transformed_reference_pose_flag": task["task_type"] != "redocking_native_receptor",
            "pocket_alignment_rmsd_for_reference": transform["pocket_alignment_rmsd"],
            "rmsd_status": status,
            "rmsd_warnings_json": json.dumps([warning] if warning else []),
        })
    frame = pd.DataFrame(rows, columns=RMSD_COLUMNS)
    write_table(paths["processed"] / "pose_rmsd.parquet", frame)
    write_table(paths["processed"] / "pose_rmsd.csv", frame)
    return frame
