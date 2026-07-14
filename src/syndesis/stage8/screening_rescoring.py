from __future__ import annotations

import json
import subprocess
import time
from pathlib import Path

import pandas as pd

from syndesis.common.io import project_root
from syndesis.stage4.gnina_runner import _prepare_ligand_for_gnina, build_gnina_command
from syndesis.stage4.score_parser import parse_gnina_output


def _rel(path: str | Path) -> str:
    p = Path(path)
    try:
        return str(p.relative_to(project_root()))
    except ValueError:
        return str(p)


def rescore_screening_poses(poses: pd.DataFrame, triage: pd.DataFrame, tasks: pd.DataFrame, config: dict, paths: dict) -> pd.DataFrame:
    receptor = tasks.set_index("screening_task_id")["receptor_file"].to_dict()
    keep = poses.merge(triage[["screening_pose_id", "pass_primary_docking_triage", "docking_score_percentile_within_receptor"]], on="screening_pose_id")
    keep = keep[keep["pass_primary_docking_triage"]]
    rows = []
    for pose in keep.to_dict("records"):
        ligand = _prepare_ligand_for_gnina(pose["pose_file"], {"processed": paths["processed"]}, config["gnina"].get("obabel") or config.get("prep", {}).get("obabel_path"))
        cmd = build_gnina_command(config["gnina"], ["--score_only", "-r", _rel(receptor[pose["screening_task_id"]]), "-l", _rel(ligand)])
        start = time.time()
        proc = subprocess.run(cmd, check=False, capture_output=True, text=True, timeout=int(config["gnina"]["timeout_seconds"]), cwd=project_root())
        parsed = parse_gnina_output(proc.stdout)
        rows.append({"screening_pose_id": pose["screening_pose_id"], "molecule_id": pose["molecule_id"], "target_receptor_id": pose["target_receptor_id"], "receptor_state": pose["receptor_state"], "docking_score": pose["docking_score"], "docking_score_percentile_within_receptor": pose["docking_score_percentile_within_receptor"], "gnina_empirical_affinity": parsed["gnina_empirical_affinity"], "cnnscore": parsed["cnnscore"], "cnnaffinity": parsed["cnnaffinity"], "cnn_vs": parsed["cnn_vs"], "gnina_model": config["gnina"]["model_name"], "gnina_version": "gnina", "rescoring_status": "success" if proc.returncode == 0 else "failed", "warnings_json": json.dumps([] if proc.returncode == 0 else [proc.stderr[-500:]])})
    out = pd.DataFrame(rows)
    out.to_parquet(paths["processed"] / "screening_gnina_scores.parquet", index=False)
    out.to_csv(paths["processed"] / "screening_gnina_scores.csv", index=False)
    return out
