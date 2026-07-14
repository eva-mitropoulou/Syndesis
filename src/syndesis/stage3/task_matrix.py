from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from Bio.PDB import PDBParser, Superimposer

from syndesis.common.io import ensure_dir, write_table
from syndesis.stage3.schemas import TASK_COLUMNS, TRANSFORM_COLUMNS
from syndesis.stage2.pocket_mapping import resolve_uniprot_residue, residue_by_auth_seq


def parse_json_list(value: str | None) -> list[float]:
    if not isinstance(value, str):
        return [None, None, None]
    return list(json.loads(value))


def write_transform_matrix(path: Path) -> Path:
    ensure_dir(path.parent)
    np.savetxt(path, np.eye(4), fmt="%.8f")
    return path


def ca_atoms(path: str | Path, residues: list[int]) -> dict[int, object]:
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure(Path(path).stem, str(path))
    chain = next(structure.get_chains())
    by_seq = residue_by_auth_seq(chain)
    atoms = {}
    for residue_number in residues:
        residue, _, _ = resolve_uniprot_residue(by_seq, int(residue_number))
        if residue is None:
            continue
        for atom in residue.get_atoms():
            if atom.get_name().strip() == "CA":
                atoms[int(residue_number)] = atom
    return atoms


def pocket_transform(native_receptor: str | Path, target_receptor: str | Path, residues: list[int]) -> tuple[np.ndarray, np.ndarray, float | None, int, list[str]]:
    native_atoms = ca_atoms(native_receptor, residues)
    target_atoms = ca_atoms(target_receptor, residues)
    shared = sorted(set(native_atoms) & set(target_atoms))
    warnings: list[str] = []
    if len(shared) < 3:
        warnings.append("Fewer than 3 shared pocket C-alpha atoms; identity transform used.")
        return np.eye(3), np.zeros(3), None, len(shared), warnings
    sup = Superimposer()
    sup.set_atoms([target_atoms[i] for i in shared], [native_atoms[i] for i in shared])
    rot, tran = sup.rotran
    return np.asarray(rot), np.asarray(tran), float(sup.rms), len(shared), warnings


def transform_pdb_coords(input_pdb: str | Path, output_pdb: Path, rot: np.ndarray, tran: np.ndarray) -> None:
    ensure_dir(output_pdb.parent)
    lines = []
    with Path(input_pdb).open("r", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith(("ATOM", "HETATM")):
                coord = np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])
                new = np.dot(coord, rot) + tran
                line = f"{line[:30]}{new[0]:8.3f}{new[1]:8.3f}{new[2]:8.3f}{line[54:]}"
            lines.append(line)
    output_pdb.write_text("".join(lines), encoding="utf-8")


def write_matrix(path: Path, rot: np.ndarray, tran: np.ndarray) -> Path:
    ensure_dir(path.parent)
    matrix = np.eye(4)
    matrix[:3, :3] = rot
    matrix[3, :3] = tran
    np.savetxt(path, matrix, fmt="%.8f")
    return path


def build_reference_transforms(ensemble: pd.DataFrame, ligand_prep: pd.DataFrame, paths: dict[str, Path]) -> pd.DataFrame:
    rows = []
    for _, ligand in ligand_prep.iterrows():
        native = ensemble[ensemble["complex_id"] == ligand["complex_id"]].iloc[0]
        for _, target in ensemble.iterrows():
            matrix_path = paths["reference_transforms"] / f"{ligand['ligand_id']}__{native['receptor_id']}__{target['receptor_id']}.txt"
            transformed = paths["reference_transforms"] / f"{ligand['ligand_id']}__{target['receptor_id']}_reference.pdbqt"
            residues = [745, 762, 790, 793, 797, 855, 856, 857]
            if native["receptor_id"] == target["receptor_id"]:
                rot, tran, align_rmsd, align_count, warnings = np.eye(3), np.zeros(3), 0.0, len(residues), []
            else:
                rot, tran, align_rmsd, align_count, warnings = pocket_transform(native["receptor_file_path"], target["receptor_file_path"], residues)
            transform_pdb_coords(ligand["docking_format_file"], transformed, rot, tran)
            rows.append({
                "ligand_id": ligand["ligand_id"],
                "native_receptor_id": native["receptor_id"],
                "target_receptor_id": target["receptor_id"],
                "transform_matrix_file": str(write_matrix(matrix_path, rot, tran)),
                "pocket_alignment_rmsd": align_rmsd,
                "alignment_residue_count": align_count,
                "transformed_reference_pose_file": str(transformed),
                "transform_status": "identity" if native["receptor_id"] == target["receptor_id"] else "pocket_aligned",
                "warnings_json": json.dumps(warnings),
            })
    frame = pd.DataFrame(rows, columns=TRANSFORM_COLUMNS)
    write_table(paths["processed"] / "reference_pose_transforms.parquet", frame)
    write_table(paths["processed"] / "reference_pose_transforms.csv", frame)
    return frame


def task_type(native_state: str, target_state: str, native_receptor_id: str, target_receptor_id: str) -> tuple[str, bool]:
    if native_receptor_id == target_receptor_id:
        return "redocking_native_receptor", True
    state_match = native_state == target_state and native_state != "unknown_state"
    return ("crossdocking_same_state" if state_match else "crossdocking_other_state"), state_match


def build_task_matrix(
    ensemble: pd.DataFrame,
    receptor_prep: pd.DataFrame,
    ligand_prep: pd.DataFrame,
    transforms: pd.DataFrame,
    config: dict[str, Any],
    paths: dict[str, Path],
) -> pd.DataFrame:
    receptor_lookup = receptor_prep.set_index("receptor_id")
    transform_lookup = transforms.set_index(["ligand_id", "target_receptor_id"])
    rows = []
    engine = config["docking"]["primary_engine"]
    seed = int(config["docking"]["seed_list"][0])
    for _, ligand in ligand_prep.iterrows():
        native = ensemble[ensemble["complex_id"] == ligand["complex_id"]].iloc[0]
        for _, target in ensemble.iterrows():
            ttype, state_match = task_type(native["state_stratum"], target["state_stratum"], native["receptor_id"], target["receptor_id"])
            center = parse_json_list(target["suggested_docking_box_center"])
            size = parse_json_list(target["suggested_docking_box_size"])
            transform = transform_lookup.loc[(ligand["ligand_id"], target["receptor_id"])]
            task_id = f"{ligand['ligand_id']}__{target['receptor_id']}__{engine}__seed{seed}"
            rows.append({
                "docking_task_id": task_id,
                "task_type": ttype,
                "ligand_id": ligand["ligand_id"],
                "ligand_native_complex_id": ligand["complex_id"],
                "native_receptor_id": native["receptor_id"],
                "target_receptor_id": target["receptor_id"],
                "target_receptor_state": target["state_stratum"],
                "native_receptor_state": native["state_stratum"],
                "state_match_flag": state_match,
                "ligand_prepared_file": ligand["docking_format_file"],
                "receptor_prepared_file": receptor_lookup.loc[target["receptor_id"], "docking_format_file"],
                "docking_box_center_x": center[0],
                "docking_box_center_y": center[1],
                "docking_box_center_z": center[2],
                "docking_box_size_x": size[0],
                "docking_box_size_y": size[1],
                "docking_box_size_z": size[2],
                "docking_engine": engine,
                "engine_version": "unavailable",
                "exhaustiveness": config["docking"]["exhaustiveness"],
                "num_modes": config["docking"]["num_modes"],
                "seed": seed,
                "replicate_id": 1,
                "expected_native_pose_reference_file": transform["transformed_reference_pose_file"],
                "reference_pose_transform_matrix_path": transform["transform_matrix_file"],
                "task_status": "pending",
                "skip_reason": "",
            })
    frame = pd.DataFrame(rows, columns=TASK_COLUMNS)
    write_table(paths["processed"] / "docking_task_matrix.parquet", frame)
    write_table(paths["processed"] / "docking_task_matrix.csv", frame)
    return frame
