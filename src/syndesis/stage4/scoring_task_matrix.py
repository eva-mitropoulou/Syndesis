from __future__ import annotations

from pathlib import Path
from typing import Any

import pandas as pd

from syndesis.common.io import write_table
from syndesis.stage4.schemas import RESCORING_TASK_COLUMNS


def _exists(value: object) -> bool:
    return isinstance(value, str) and bool(value) and Path(value).exists()


def build_rescoring_tasks(inputs: dict[str, pd.DataFrame], config: dict[str, Any], paths: dict[str, Path]) -> pd.DataFrame:
    poses = inputs["poses"].rename(columns={"pose_rank": "original_pose_rank", "docking_score": "original_docking_score"})
    tasks = inputs["tasks"][["docking_task_id", "task_type", "receptor_prepared_file", "target_receptor_state"]]
    labels = inputs["labels"][
        [
            "pose_id", "stage3_pose_label", "rmsd_symmetry_corrected", "sanity_status",
            "invalid_pose_flag",
        ]
    ]
    merged = poses.merge(tasks, on="docking_task_id", how="left").merge(labels, on="pose_id", how="left")
    top_n = int(config["task_matrix"].get("top_n_per_docking_task", 20))
    include_invalid = bool(config["task_matrix"].get("include_invalid_poses", False))
    merged = merged[merged["original_pose_rank"] <= top_n].copy()

    rows: list[dict[str, Any]] = []
    for idx, row in merged.iterrows():
        warnings: list[str] = []
        status = "ready"
        skip_reason = ""
        if not _exists(row.get("pose_file")):
            status = "skipped"
            skip_reason = "missing_pose_file"
        elif not _exists(row.get("receptor_prepared_file")):
            status = "skipped"
            skip_reason = "missing_receptor_file"
        elif bool(row.get("invalid_pose_flag")) and not include_invalid:
            status = "skipped"
            skip_reason = "invalid_pose_excluded_by_config"
        if pd.isna(row.get("stage3_pose_label")):
            warnings.append("missing_stage3_label")
        rows.append(
            {
                "rescoring_task_id": f"rescore__{row['pose_id']}",
                "pose_id": row["pose_id"],
                "docking_task_id": row["docking_task_id"],
                "ligand_id": row["ligand_id"],
                "target_receptor_id": row["target_receptor_id"],
                "task_type": row.get("task_type"),
                "docking_engine": row["docking_engine"],
                "original_pose_rank": int(row["original_pose_rank"]),
                "original_docking_score": row["original_docking_score"],
                "pose_file": row["pose_file"],
                "receptor_file": row.get("receptor_prepared_file"),
                "receptor_state": row.get("target_receptor_state"),
                "native_like_label_stage3": row.get("stage3_pose_label"),
                "rmsd_symmetry_corrected": row.get("rmsd_symmetry_corrected"),
                "sanity_status": row.get("sanity_status"),
                "rescoring_engine": "gnina",
                "rescoring_mode": "score_only",
                "model_name": config["gnina"].get("model_name", "default"),
                "model_version": config["gnina"].get("model_version"),
                "use_gpu": bool(config["gnina"].get("use_gpu", True)),
                "batch_id": idx // max(int(config["gnina"].get("batch_size", 1)), 1),
                "task_status": status,
                "skip_reason": skip_reason,
            }
        )
    frame = pd.DataFrame(rows, columns=RESCORING_TASK_COLUMNS)
    write_table(paths["processed"] / "rescoring_task_matrix.parquet", frame)
    write_table(paths["processed"] / "rescoring_task_matrix.csv", frame)
    return frame
