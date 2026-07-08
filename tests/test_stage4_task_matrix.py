from __future__ import annotations

from pathlib import Path

import pandas as pd

from egfr_dockingforge.stage4.scoring_task_matrix import build_rescoring_tasks


def test_rescoring_task_matrix_skips_missing_pose(tmp_path: Path) -> None:
    receptor = tmp_path / "receptor.pdb"
    receptor.write_text("END\n", encoding="utf-8")
    inputs = {
        "poses": pd.DataFrame(
            [
                {
                    "pose_id": "pose1",
                    "docking_task_id": "task1",
                    "ligand_id": "lig1",
                    "target_receptor_id": "rec1",
                    "docking_engine": "unidock",
                    "pose_rank": 1,
                    "docking_score": -7.0,
                    "pose_file": str(tmp_path / "missing.pdbqt"),
                }
            ]
        ),
        "tasks": pd.DataFrame(
            [{"docking_task_id": "task1", "task_type": "redocking_native_receptor", "receptor_prepared_file": str(receptor), "target_receptor_state": "active_like"}]
        ),
        "labels": pd.DataFrame(
            [{"pose_id": "pose1", "stage3_pose_label": "sampling_failure", "rmsd_symmetry_corrected": 8.0, "sanity_status": "pass", "invalid_pose_flag": False}]
        ),
    }
    config = {"task_matrix": {"top_n_per_docking_task": 20, "include_invalid_poses": False}, "gnina": {"model_name": "default", "model_version": None, "use_gpu": True, "batch_size": 1}}
    paths = {"processed": tmp_path}
    tasks = build_rescoring_tasks(inputs, config, paths)
    assert tasks.loc[0, "task_status"] == "skipped"
    assert tasks.loc[0, "skip_reason"] == "missing_pose_file"

