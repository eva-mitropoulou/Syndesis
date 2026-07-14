from __future__ import annotations

from syndesis.stage4.schemas import GNINA_SCORE_COLUMNS, POSE_SCORE_COLUMNS, RESCORING_TASK_COLUMNS, TASK_METRIC_COLUMNS


def test_stage4_required_schema_columns_exist() -> None:
    assert {"rescoring_task_id", "pose_id", "receptor_file", "task_status"}.issubset(RESCORING_TASK_COLUMNS)
    assert {"cnnscore", "cnnaffinity", "gnina_empirical_affinity"}.issubset(GNINA_SCORE_COLUMNS)
    assert {"interaction_recovery_status", "final_pose_label_status"}.issubset(POSE_SCORE_COLUMNS)
    assert {"top1_rmsd_cnnscore", "spearman_score_vs_rmsd_cnnscore"}.issubset(TASK_METRIC_COLUMNS)

