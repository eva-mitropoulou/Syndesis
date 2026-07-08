from __future__ import annotations

import json
import shutil
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from Bio.PDB import PDBIO, Superimposer

from egfr_dockingforge.common.io import ensure_dir
from egfr_dockingforge.stage2.pocket_mapping import first_chain, parse_receptor, residue_by_auth_seq
from egfr_dockingforge.stage2.schemas import ALIGNMENT_COLUMNS


def ca_atoms_for_residues(path: str | Path, residues: list[int]) -> dict[int, Any]:
    structure = parse_receptor(path)
    chain = first_chain(structure)
    by_seq = residue_by_auth_seq(chain)
    atoms = {}
    for residue_number in residues:
        residue = by_seq.get(int(residue_number))
        if residue is None:
            continue
        for atom in residue.get_atoms():
            if atom.get_name().strip() == "CA":
                atoms[int(residue_number)] = atom
    return atoms


def rmsd(coords_a: np.ndarray, coords_b: np.ndarray) -> float | None:
    if len(coords_a) == 0 or len(coords_a) != len(coords_b):
        return None
    return float(np.sqrt(((coords_a - coords_b) ** 2).sum(axis=1).mean()))


def receptor_pair_rmsd(path_a: str | Path, path_b: str | Path, residues: list[int]) -> tuple[float | None, int, int]:
    atoms_a = ca_atoms_for_residues(path_a, residues)
    atoms_b = ca_atoms_for_residues(path_b, residues)
    shared = sorted(set(atoms_a) & set(atoms_b))
    missing = len(set(residues) - set(shared))
    if len(shared) < 3:
        return None, len(shared), missing
    sup = Superimposer()
    sup.set_atoms([atoms_a[i] for i in shared], [atoms_b[i] for i in shared])
    return float(sup.rms), len(shared), missing


def superimposer_for_receptor(
    reference_path: str | Path, mobile_path: str | Path, residues: list[int]
) -> Superimposer | None:
    """Return a fitted Superimposer mapping the mobile receptor onto the reference.

    Uses shared pocket CA atoms (same set used for the RMSD metric). Returns None when
    fewer than three shared CA atoms are available, in which case no superposition can be
    defined and the mobile coordinates must be left untransformed.
    """
    ref_atoms = ca_atoms_for_residues(reference_path, residues)
    mobile_atoms = ca_atoms_for_residues(mobile_path, residues)
    shared = sorted(set(ref_atoms) & set(mobile_atoms))
    if len(shared) < 3:
        return None
    sup = Superimposer()
    # set_atoms(fixed, moving): rotate/translate the moving (mobile) atoms onto the fixed
    # (reference) frame; sup.rotran can then be applied to the whole mobile structure.
    sup.set_atoms([ref_atoms[i] for i in shared], [mobile_atoms[i] for i in shared])
    return sup


def copy_aligned_receptors(features: pd.DataFrame, out_dir: Path, config: dict[str, Any]) -> pd.DataFrame:
    """Superpose each receptor onto the ensemble reference and write the aligned PDB.

    Previously this copied the ORIGINAL receptor coordinates and labelled them "aligned",
    discarding the computed Superimposer transform. We now apply the pocket-CA superposition
    (rot, tran) to every receptor's atoms before writing, so the output files are genuinely
    in the reference frame. When a superposition cannot be defined (fewer than three shared
    pocket CA atoms) we fall back to an untransformed copy and record a warning so downstream
    consumers do not mistake it for a superposed structure.
    """
    ensure_dir(out_dir)
    residues = [int(x) for x in config["pocket"]["alignment_residues_uniprot"]]
    rows = []
    if features.empty:
        return pd.DataFrame(rows, columns=["receptor_id", "aligned_receptor_file_path", "warnings_json"])
    reference = features.sort_values(["quality_tier", "resolution_angstrom", "receptor_id"]).iloc[0]
    reference_path = reference["receptor_file_path"]
    io = PDBIO()
    for _, row in features.iterrows():
        target = out_dir / f"{row['receptor_id']}.pdb"
        warnings: list[str] = []
        if row["receptor_id"] == reference["receptor_id"]:
            # The reference defines the frame; copy it verbatim.
            shutil.copyfile(row["receptor_file_path"], target)
        else:
            sup = superimposer_for_receptor(reference_path, row["receptor_file_path"], residues)
            if sup is None:
                shutil.copyfile(row["receptor_file_path"], target)
                warnings.append("receptor_not_superposed_fewer_than_3_shared_pocket_ca_atoms")
            else:
                structure = parse_receptor(row["receptor_file_path"])
                rot, tran = sup.rotran
                structure.transform(rot, tran)
                io.set_structure(structure)
                io.save(str(target))
        rows.append(
            {
                "receptor_id": row["receptor_id"],
                "aligned_receptor_file_path": str(target),
                "warnings_json": json.dumps(warnings),
            }
        )
    return pd.DataFrame(rows, columns=["receptor_id", "aligned_receptor_file_path", "warnings_json"])


def alignment_metrics(features: pd.DataFrame, config: dict[str, Any]) -> pd.DataFrame:
    residues = [int(x) for x in config["pocket"]["alignment_residues_uniprot"]]
    if features.empty:
        return pd.DataFrame(columns=ALIGNMENT_COLUMNS)
    reference = features.sort_values(["quality_tier", "resolution_angstrom", "receptor_id"]).iloc[0]
    rows = []
    for _, row in features.iterrows():
        pocket_rmsd, shared, missing = receptor_pair_rmsd(reference["receptor_file_path"], row["receptor_file_path"], residues)
        hinge_rmsd, _, _ = receptor_pair_rmsd(reference["receptor_file_path"], row["receptor_file_path"], [790, 793, 797])
        dfg_rmsd, _, _ = receptor_pair_rmsd(reference["receptor_file_path"], row["receptor_file_path"], [855, 856, 857])
        rows.append(
            {
                "receptor_id": row["receptor_id"],
                "reference_receptor_id": reference["receptor_id"],
                "pocket_ca_rmsd": pocket_rmsd,
                "pocket_backbone_rmsd": pocket_rmsd,
                "pocket_sidechain_rmsd": None,
                "hinge_region_rmsd": hinge_rmsd,
                "dfg_region_rmsd": dfg_rmsd,
                "activation_loop_rmsd_if_available": None,
                "alignment_residue_count": shared,
                "missing_alignment_residue_count": missing,
            }
        )
    return pd.DataFrame(rows, columns=ALIGNMENT_COLUMNS)


def update_features_with_alignment(features: pd.DataFrame, metrics: pd.DataFrame) -> pd.DataFrame:
    if metrics.empty:
        return features
    mapped = metrics.set_index("receptor_id")
    features = features.copy()
    for column, metric in [
        ("pocket_ca_rmsd_to_reference", "pocket_ca_rmsd"),
        ("hinge_region_rmsd", "hinge_region_rmsd"),
        ("dfg_region_rmsd", "dfg_region_rmsd"),
    ]:
        features[column] = features["receptor_id"].map(mapped[metric])
    return features
