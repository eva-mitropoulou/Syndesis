from __future__ import annotations

from pathlib import Path

import pandas as pd

from syndesis.stage4.score_diagnostics import failure_table, rescoring_task_metrics


def test_ranking_directions_are_respected(tmp_path: Path) -> None:
    table = pd.DataFrame(
        [
            {"pose_id": "p1", "docking_task_id": "t", "task_type": "redocking_native_receptor", "ligand_id": "l", "target_receptor_id": "r", "docking_engine": "unidock", "original_pose_rank": 1, "original_docking_score": -8.0, "cnnscore": 0.1, "cnnaffinity": 1.0, "gnina_empirical_affinity": -6.0, "rmsd_symmetry_corrected": 8.0, "strict_native_like_flag": False, "relaxed_native_like_flag": False, "invalid_pose_flag": False, "rescoring_status": "success"},
            {"pose_id": "p2", "docking_task_id": "t", "task_type": "redocking_native_receptor", "ligand_id": "l", "target_receptor_id": "r", "docking_engine": "unidock", "original_pose_rank": 2, "original_docking_score": -6.0, "cnnscore": 0.9, "cnnaffinity": 4.0, "gnina_empirical_affinity": -7.0, "rmsd_symmetry_corrected": 1.0, "strict_native_like_flag": True, "relaxed_native_like_flag": True, "invalid_pose_flag": False, "rescoring_status": "success"},
        ]
    )
    config = {"ranking": {"strict_rmsd_threshold_angstrom": 2.0, "relaxed_rmsd_threshold_angstrom": 3.0, "directions": {"original_docking_score": "lower", "cnnscore": "higher", "cnnaffinity": "higher", "gnina_empirical_affinity": "lower"}}}
    metrics = rescoring_task_metrics(table, config, {"processed": tmp_path})
    assert metrics.loc[0, "top1_rmsd_original_score"] == 8.0
    assert metrics.loc[0, "top1_rmsd_cnnscore"] == 1.0
    assert bool(metrics.loc[0, "top1_success_cnnscore_strict"]) is True


def test_failure_table_detects_wrong_cnn_preference(tmp_path: Path) -> None:
    table = pd.DataFrame(
        [
            {"pose_id": "p1", "docking_task_id": "t", "task_type": "redocking_native_receptor", "ligand_id": "l", "target_receptor_id": "r", "original_docking_score": -8.0, "cnnscore": 0.9, "cnnaffinity": 3.0, "rmsd_symmetry_corrected": 8.0, "sanity_status": "pass", "stage3_pose_label": "sampling_failure", "relaxed_native_like_flag": False, "rescoring_status": "success"},
            {"pose_id": "p2", "docking_task_id": "t", "task_type": "redocking_native_receptor", "ligand_id": "l", "target_receptor_id": "r", "original_docking_score": -6.0, "cnnscore": 0.1, "cnnaffinity": 1.0, "rmsd_symmetry_corrected": 1.0, "sanity_status": "pass", "stage3_pose_label": "strict_native_like", "relaxed_native_like_flag": True, "rescoring_status": "success"},
        ]
    )
    failures = failure_table(table, {"processed": tmp_path})
    assert "cnnscore_prefers_wrong_pose" in set(failures["failure_type"])

