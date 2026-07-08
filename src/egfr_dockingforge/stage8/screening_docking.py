from __future__ import annotations

import hashlib
import json
import subprocess
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import pandas as pd

from egfr_dockingforge.common.io import write_table
from egfr_dockingforge.stage3.unidock_runner import split_pdbqt_models


def _now() -> str:
    return datetime.now(timezone.utc).isoformat()


def run_candidate_docking(tasks: pd.DataFrame, config: dict[str, Any], paths: dict[str, Path]) -> tuple[pd.DataFrame, pd.DataFrame]:
    runs = []
    poses = []
    cfg_hash = hashlib.sha256(json.dumps(config, sort_keys=True, default=str).encode()).hexdigest()[:16]
    for task in tasks.to_dict("records"):
        run_id = f"run__{task['screening_task_id']}"
        raw = paths["docking_outputs_raw"] / f"{run_id}.pdbqt"
        log = paths["docking_outputs_logs"] / f"{run_id}.log"
        cmd = [config["docking"]["executable"], "--receptor", task["receptor_file"], "--ligand", task["ligand_file"], "--center_x", str(task["docking_box_center_x"]), "--center_y", str(task["docking_box_center_y"]), "--center_z", str(task["docking_box_center_z"]), "--size_x", str(task["docking_box_size_x"]), "--size_y", str(task["docking_box_size_y"]), "--size_z", str(task["docking_box_size_z"]), "--exhaustiveness", str(task["exhaustiveness"]), "--num_modes", str(task["num_modes"]), "--seed", str(task["seed"]), "--cpu", "4", "--out", str(raw)]
        start = time.time()
        start_iso = _now()
        completed = subprocess.run(cmd, check=False, capture_output=True, text=True, timeout=int(config["docking"]["timeout_seconds"]))
        log.write_text(completed.stdout + "\n" + completed.stderr, encoding="utf-8")
        status = "success" if completed.returncode == 0 and raw.exists() else "failed"
        runs.append({"screening_run_id": run_id, "screening_task_id": task["screening_task_id"], "docking_engine": task["docking_engine"], "engine_version": task["engine_version"], "command_line": " ".join(cmd), "start_time": start_iso, "end_time": _now(), "runtime_seconds": round(time.time() - start, 3), "exit_code": completed.returncode, "status": status, "error_message": "" if status == "success" else completed.stderr[-1000:], "output_pose_file": str(raw) if raw.exists() else "", "output_log_file": str(log), "config_hash": cfg_hash})
        if status == "success":
            for pose in split_pdbqt_models(raw, paths["docking_outputs_poses"], run_id):
                sid = f"{run_id}__pose{pose['pose_rank']:02d}"
                poses.append({"screening_pose_id": sid, "screening_run_id": run_id, "screening_task_id": task["screening_task_id"], "prepared_ligand_id": task["prepared_ligand_id"], "molecule_id": task["molecule_id"], "source": task["source"], "target_receptor_id": task["target_receptor_id"], "receptor_state": task["receptor_state"], "docking_engine": task["docking_engine"], "protonation_state_id": "stage7_rdkit_default", "tautomer_state_id": "canonical_tautomer", "conformer_id": "conf_0", "seed": task["seed"], "pose_rank": pose["pose_rank"], "docking_score": pose["docking_score"], "pose_file": pose["pose_file"], "parse_status": "complete", "parse_warnings_json": json.dumps([])})
    run_df = pd.DataFrame(runs)
    pose_df = pd.DataFrame(poses)
    write_table(paths["processed"] / "screening_docking_runs.parquet", run_df)
    write_table(paths["processed"] / "screening_docking_runs.csv", run_df)
    write_table(paths["processed"] / "screening_docked_poses.parquet", pose_df)
    write_table(paths["processed"] / "screening_docked_poses.csv", pose_df)
    return run_df, pose_df
