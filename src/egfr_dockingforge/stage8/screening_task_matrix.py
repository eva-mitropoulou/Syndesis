from __future__ import annotations

import ast
import json
import subprocess
from pathlib import Path
from typing import Any

import pandas as pd

from egfr_dockingforge.common.io import write_table


def _vec(value: str | list[float]) -> list[float]:
    return [float(x) for x in (ast.literal_eval(value) if isinstance(value, str) else value)]


def prepare_manifest(stage7_input: pd.DataFrame, master: pd.DataFrame, config: dict[str, Any], paths: dict[str, Path]) -> pd.DataFrame:
    merged = stage7_input.merge(master[["molecule_id", "subsource", "hard_scope_pass", "include_in_screening_library", "risk_flags_json"]], on="molecule_id", how="left")
    controls = merged[merged["screening_role"].isin(["known_activity_reference", "native_pose_reference"])].head(int(config["screening"]["max_known_controls"]))
    candidates = merged[merged["include_in_screening_library"].fillna(False).astype(bool)].head(int(config["screening"]["max_candidates"]))
    chosen = pd.concat([controls, candidates], ignore_index=True).drop_duplicates("prepared_ligand_id")
    rows = []
    for row in chosen.to_dict("records"):
        ligand_file = row.get("sdf_path")
        include = bool(ligand_file and Path(ligand_file).exists() and row.get("source") and row.get("screening_role") and row.get("hard_scope_pass", True))
        rows.append({"prepared_ligand_id": row["prepared_ligand_id"], "molecule_id": row["molecule_id"], "source": row["source"], "subsource": row.get("subsource", ""), "screening_role": row["screening_role"], "screening_subset": row["screening_subset"], "standard_smiles": row["standard_smiles"], "prepared_smiles": row["prepared_smiles"], "ligand_file": ligand_file, "novelty_bucket": row["novelty_bucket"], "closest_known_molecule_id": row["closest_known_molecule_id"], "tanimoto_to_closest_known": row["tanimoto_to_closest_known"], "medchem_flags_json": row.get("risk_flags_json", "[]"), "include_in_screening": include, "exclusion_reason": "" if include else "missing_ligand_or_hard_scope_failed", "warnings_json": json.dumps([])})
    manifest = pd.DataFrame(rows)
    write_table(paths["processed"] / "screening_candidate_manifest.parquet", manifest)
    write_table(paths["processed"] / "screening_candidate_manifest.csv", manifest)
    return manifest


def _prepare_pdbqt(sdf_path: str, prepared_id: str, config: dict[str, Any], paths: dict[str, Path]) -> str:
    target = paths["prepared_ligands"] / f"{prepared_id}.pdbqt"
    if target.exists() and target.stat().st_mtime >= Path(sdf_path).stat().st_mtime:
        return str(target)
    cmd = [config["docking"]["obabel"], "-isdf", sdf_path, "-opdbqt", "-O", str(target), "--partialcharge", "gasteiger"]
    completed = subprocess.run(cmd, check=False, capture_output=True, text=True, timeout=120)
    if completed.returncode != 0 or not target.exists():
        raise RuntimeError(f"OpenBabel PDBQT conversion failed for {sdf_path}: {completed.stderr[-500:]}")
    return str(target)


def build_screening_task_matrix(manifest: pd.DataFrame, receptors: pd.DataFrame, config: dict[str, Any], paths: dict[str, Path]) -> pd.DataFrame:
    selected = receptors[receptors["selected_flag"].fillna(False).astype(bool)].head(int(config["screening"]["receptor_limit"]))
    rows = []
    for ligand in manifest[manifest["include_in_screening"]].to_dict("records"):
        pdbqt = _prepare_pdbqt(ligand["ligand_file"], ligand["prepared_ligand_id"], config, paths)
        for receptor in selected.to_dict("records"):
            center = _vec(receptor["suggested_docking_box_center"])
            size = _vec(receptor["suggested_docking_box_size"])
            task_id = f"screen__{ligand['prepared_ligand_id']}__{receptor['receptor_id']}"
            rows.append({"screening_task_id": task_id, "prepared_ligand_id": ligand["prepared_ligand_id"], "molecule_id": ligand["molecule_id"], "source": ligand["source"], "screening_role": ligand["screening_role"], "screening_subset": ligand["screening_subset"], "target_receptor_id": receptor["receptor_id"], "receptor_state": receptor["state_stratum"], "receptor_cluster_id": receptor["cluster_id"], "receptor_file": receptor["receptor_file_path"], "ligand_file": pdbqt, "docking_box_center_x": center[0], "docking_box_center_y": center[1], "docking_box_center_z": center[2], "docking_box_size_x": size[0], "docking_box_size_y": size[1], "docking_box_size_z": size[2], "docking_engine": config["docking"]["engine"], "engine_version": "available", "num_modes": int(config["docking"]["num_modes"]), "exhaustiveness": int(config["docking"]["exhaustiveness"]), "seed": int(config["docking"]["seed"]), "batch_id": "smoke_batch_0", "task_status": "ready", "skip_reason": ""})
    tasks = pd.DataFrame(rows)
    write_table(paths["processed"] / "screening_task_matrix.parquet", tasks)
    write_table(paths["processed"] / "screening_task_matrix.csv", tasks)
    return tasks
