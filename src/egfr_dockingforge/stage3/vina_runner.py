from __future__ import annotations

import json
import hashlib
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import pandas as pd

from egfr_dockingforge.common.io import ensure_dir, write_table
from egfr_dockingforge.stage3.docking_engines import check_engine
from egfr_dockingforge.stage3.schemas import DOCKED_POSE_COLUMNS, DOCKING_RUN_COLUMNS


def now() -> str:
    return datetime.now(timezone.utc).isoformat()


def config_hash(config: dict[str, Any]) -> str:
    payload = json.dumps(config, sort_keys=True, default=str)
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def run_vina_tasks(task_matrix: pd.DataFrame, config: dict[str, Any], paths: dict[str, Path]) -> tuple[pd.DataFrame, pd.DataFrame]:
    engine = config["docking"]["primary_engine"]
    executable = config["docking"].get("engine_executable_path", {}).get(engine, engine)
    availability = check_engine(engine, executable)
    cfg_hash = config_hash(config)
    rows = []
    pose_rows = []
    ensure_dir(paths["docking_outputs_logs"])
    for _, task in task_matrix.iterrows():
        run_id = f"run__{task['docking_task_id']}"
        log_path = paths["docking_outputs_logs"] / f"{run_id}.log"
        command = (
            f"{executable} --receptor {task['receptor_prepared_file']} --ligand {task['ligand_prepared_file']} "
            f"--center_x {task['docking_box_center_x']} --center_y {task['docking_box_center_y']} --center_z {task['docking_box_center_z']} "
            f"--size_x {task['docking_box_size_x']} --size_y {task['docking_box_size_y']} --size_z {task['docking_box_size_z']} "
            f"--exhaustiveness {task['exhaustiveness']} --num_modes {task['num_modes']} --seed {task['seed']}"
        )
        status = "engine_unavailable" if not availability.available else "not_run_manual_execution_required"
        error = None if availability.available else f"Docking engine executable not found: {executable}"
        log_path.write_text(json.dumps({"status": status, "error": error, "command_line": command}, indent=2) + "\n", encoding="utf-8")
        start = now()
        rows.append({
            "docking_run_id": run_id,
            "docking_task_id": task["docking_task_id"],
            "docking_engine": engine,
            "engine_version": availability.version,
            "command_line": command,
            "container_image": config["docking"].get("container_image", {}).get(engine),
            "start_time": start,
            "end_time": start,
            "runtime_seconds": 0.0,
            "exit_code": 127 if not availability.available else None,
            "status": status,
            "error_message": error,
            "output_pose_file": None,
            "output_log_file": str(log_path),
            "config_hash": cfg_hash,
        })
    runs = pd.DataFrame(rows, columns=DOCKING_RUN_COLUMNS)
    poses = pd.DataFrame(pose_rows, columns=DOCKED_POSE_COLUMNS)
    write_table(paths["processed"] / "docking_runs.parquet", runs)
    write_table(paths["processed"] / "docking_runs.csv", runs)
    write_table(paths["processed"] / "docked_poses.parquet", poses)
    write_table(paths["processed"] / "docked_poses.csv", poses)
    return runs, poses
