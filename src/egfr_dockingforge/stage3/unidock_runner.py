from __future__ import annotations

import json
import hashlib
import subprocess
import time
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
    return hashlib.sha256(json.dumps(config, sort_keys=True, default=str).encode("utf-8")).hexdigest()


def split_pdbqt_models(raw_pose_file: Path, pose_dir: Path, docking_run_id: str) -> list[dict[str, Any]]:
    if not raw_pose_file.exists():
        return []
    ensure_dir(pose_dir)
    poses: list[dict[str, Any]] = []
    current: list[str] = []
    rank = 0
    score = None
    with raw_pose_file.open("r", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("MODEL"):
                current = [line]
                rank += 1
                score = None
            elif line.startswith("ENDMDL"):
                current.append(line)
                pose_path = pose_dir / f"{docking_run_id}__pose{rank:02d}.pdbqt"
                pose_path.write_text("".join(current), encoding="utf-8")
                poses.append({"pose_rank": rank, "docking_score": score, "pose_file": str(pose_path)})
                current = []
            else:
                if current is not None:
                    current.append(line)
                if line.startswith("REMARK VINA RESULT:"):
                    parts = line.split()
                    try:
                        score = float(parts[3])
                    except Exception:
                        score = None
    return poses


def run_unidock_tasks(task_matrix: pd.DataFrame, config: dict[str, Any], paths: dict[str, Path]) -> tuple[pd.DataFrame, pd.DataFrame]:
    engine = config["docking"]["primary_engine"]
    executable = config["docking"].get("engine_executable_path", {}).get(engine, engine)
    availability = check_engine(engine, executable)
    cfg_hash = config_hash(config)
    ensure_dir(paths["docking_outputs_logs"])
    ensure_dir(paths["docking_outputs_raw"])
    ensure_dir(paths["docking_outputs_poses"])
    run_rows = []
    pose_rows = []
    for _, task in task_matrix.iterrows():
        run_id = f"run__{task['docking_task_id']}"
        raw_pose = paths["docking_outputs_raw"] / f"{run_id}.pdbqt"
        log_path = paths["docking_outputs_logs"] / f"{run_id}.log"
        command = [
            executable,
            "--receptor", str(task["receptor_prepared_file"]),
            "--ligand", str(task["ligand_prepared_file"]),
            "--center_x", str(task["docking_box_center_x"]),
            "--center_y", str(task["docking_box_center_y"]),
            "--center_z", str(task["docking_box_center_z"]),
            "--size_x", str(task["docking_box_size_x"]),
            "--size_y", str(task["docking_box_size_y"]),
            "--size_z", str(task["docking_box_size_z"]),
            "--exhaustiveness", str(task["exhaustiveness"]),
            "--num_modes", str(task["num_modes"]),
            "--seed", str(task["seed"]),
            "--cpu", "4",
            "--out", str(raw_pose),
        ]
        start = now()
        start_s = time.time()
        if not availability.available:
            status, exit_code, error = "engine_unavailable", 127, f"Docking engine executable not found: {executable}"
            log_path.write_text(json.dumps({"status": status, "error": error, "command_line": command}, indent=2), encoding="utf-8")
        else:
            completed = subprocess.run(command, check=False, capture_output=True, text=True, timeout=int(config["docking"]["timeout_seconds"]))
            log_path.write_text(completed.stdout + "\n" + completed.stderr, encoding="utf-8")
            exit_code = completed.returncode
            status = "success" if exit_code == 0 and raw_pose.exists() else "failed"
            error = None if status == "success" else completed.stderr.strip()[:1000]
        end = now()
        run_rows.append({
            "docking_run_id": run_id,
            "docking_task_id": task["docking_task_id"],
            "docking_engine": engine,
            "engine_version": availability.version,
            "command_line": " ".join(command),
            "container_image": config["docking"].get("container_image", {}).get(engine),
            "start_time": start,
            "end_time": end,
            "runtime_seconds": round(time.time() - start_s, 3),
            "exit_code": exit_code,
            "status": status,
            "error_message": error,
            "output_pose_file": str(raw_pose) if raw_pose.exists() else None,
            "output_log_file": str(log_path),
            "config_hash": cfg_hash,
        })
        if status == "success":
            for pose in split_pdbqt_models(raw_pose, paths["docking_outputs_poses"], run_id):
                pose_rank = pose["pose_rank"]
                pose_rows.append({
                    "pose_id": f"{run_id}__pose{pose_rank:02d}",
                    "docking_run_id": run_id,
                    "docking_task_id": task["docking_task_id"],
                    "ligand_id": task["ligand_id"],
                    "target_receptor_id": task["target_receptor_id"],
                    "docking_engine": engine,
                    "protonation_state_id": "native_copy",
                    "tautomer_state_id": "native_copy",
                    "conformer_id": "native_copy",
                    "seed": task["seed"],
                    "replicate_id": task["replicate_id"],
                    "pose_rank": pose_rank,
                    "docking_score": pose["docking_score"],
                    "pose_file": pose["pose_file"],
                    "raw_pose_file": str(raw_pose),
                    "parse_status": "complete",
                    "parse_warnings_json": json.dumps([]),
                })
    runs = pd.DataFrame(run_rows, columns=DOCKING_RUN_COLUMNS)
    poses = pd.DataFrame(pose_rows, columns=DOCKED_POSE_COLUMNS)
    write_table(paths["processed"] / "docking_runs.parquet", runs)
    write_table(paths["processed"] / "docking_runs.csv", runs)
    write_table(paths["processed"] / "docked_poses.parquet", poses)
    write_table(paths["processed"] / "docked_poses.csv", poses)
    return runs, poses
