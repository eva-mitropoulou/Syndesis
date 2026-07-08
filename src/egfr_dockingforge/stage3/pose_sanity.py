from __future__ import annotations

import json
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem

from egfr_dockingforge.common.io import write_table
from egfr_dockingforge.stage3.schemas import SANITY_COLUMNS

# Non-bonded heavy-atom pairs closer than this are treated as a severe internal clash.
SEVERE_CLASH_DISTANCE_ANGSTROM = 0.9
# Short (but not necessarily severe) non-bonded contacts flag a minor geometry issue.
SHORT_CONTACT_DISTANCE_ANGSTROM = 1.2


def _bonded_pairs(mol: Chem.Mol) -> set[frozenset[int]]:
    return {frozenset((bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())) for bond in mol.GetBonds()}


def _min_nonbonded_heavy_distance(mol: Chem.Mol) -> float | None:
    """Return the smallest non-bonded heavy-atom pair distance, or None if not computable."""
    conf = mol.GetConformer()
    heavy = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomicNum() > 1]
    if len(heavy) < 2:
        return None
    bonded = _bonded_pairs(mol)
    coords = {idx: np.array(conf.GetAtomPosition(idx)) for idx in heavy}
    min_dist: float | None = None
    for i, j in combinations(heavy, 2):
        if frozenset((i, j)) in bonded:
            continue
        dist = float(np.linalg.norm(coords[i] - coords[j]))
        if min_dist is None or dist < min_dist:
            min_dist = dist
    return min_dist


def sanity_for_poses(poses: pd.DataFrame, paths: dict[str, Path]) -> pd.DataFrame:
    rows = []
    for _, pose in poses.iterrows():
        pose_file = pose.get("pose_file")
        atom_loss = not isinstance(pose_file, str) or not Path(pose_file).exists()

        severe_clash = False
        ligand_geometry = False
        check_failed = False
        clash_score: float | None = None
        strain_proxy: float | None = None
        warnings: list[str] = ["PoseBusters unavailable; RDKit geometry checks only."]

        if not atom_loss:
            try:
                mol = Chem.MolFromPDBFile(pose_file, sanitize=False)
                if mol is None or mol.GetNumAtoms() == 0:
                    atom_loss = True
                elif mol.GetNumConformers() == 0:
                    ligand_geometry = True
                    warnings.append("Pose has no 3D conformer; geometry checks skipped.")
                else:
                    min_dist = _min_nonbonded_heavy_distance(mol)
                    if min_dist is not None:
                        clash_score = min_dist
                        if min_dist < SEVERE_CLASH_DISTANCE_ANGSTROM:
                            severe_clash = True
                            ligand_geometry = True
                            warnings.append(
                                f"Non-bonded heavy atoms {min_dist:.2f} A apart "
                                f"(< {SEVERE_CLASH_DISTANCE_ANGSTROM} A internal clash)."
                            )
                        elif min_dist < SHORT_CONTACT_DISTANCE_ANGSTROM:
                            ligand_geometry = True
                            warnings.append(
                                f"Closest non-bonded heavy-atom contact {min_dist:.2f} A is unusually short."
                            )
            except Exception as exc:  # noqa: BLE001 - a single bad file must not crash the batch.
                check_failed = True
                warnings.append(f"RDKit sanity check failed: {type(exc).__name__}: {exc}")

        if atom_loss or severe_clash:
            sanity_status = "failed"
        elif ligand_geometry or check_failed:
            sanity_status = "warning"
        else:
            sanity_status = "pass"

        rows.append({
            "pose_id": pose["pose_id"],
            "sanity_status": sanity_status,
            "severe_clash_flag": severe_clash,
            "ligand_geometry_flag": ligand_geometry,
            "chirality_issue_flag": False,
            "atom_loss_flag": atom_loss,
            "outside_pocket_flag": False,
            "protein_ligand_clash_score": clash_score,
            "ligand_strain_proxy": strain_proxy,
            "posebusters_available_flag": False,
            "posebusters_pass_flag": None,
            "sanity_warnings_json": json.dumps(warnings),
        })
    frame = pd.DataFrame(rows, columns=SANITY_COLUMNS)
    write_table(paths["processed"] / "pose_sanity.parquet", frame)
    write_table(paths["processed"] / "pose_sanity.csv", frame)
    return frame
