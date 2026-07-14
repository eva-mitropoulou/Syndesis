from __future__ import annotations

from syndesis.stage5 import schemas


def test_stage5_schema_exports_required_columns() -> None:
    assert "pose_id" in schemas.FINAL_POSE_LABEL_COLUMNS
    assert "final_pose_label" in schemas.FINAL_POSE_LABEL_COLUMNS
    assert "key_interaction_recall_native" in schemas.INTERACTION_RECOVERY_COLUMNS
    assert "fingerprint_sparse_json" in schemas.STAGE6_INTERACTION_FEATURE_COLUMNS

