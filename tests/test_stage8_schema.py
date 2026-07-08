from __future__ import annotations

from egfr_dockingforge.stage8 import schemas


def test_stage8_schema_has_required_columns() -> None:
    assert "prepared_ligand_id" in schemas.MANIFEST_COLUMNS
    assert "screening_task_id" in schemas.TASK_COLUMNS
    assert "pose_confidence_probability" in schemas.POSE_CONFIDENCE_COLUMNS
