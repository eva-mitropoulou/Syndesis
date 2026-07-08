from __future__ import annotations

from egfr_dockingforge.stage6 import schemas


def test_stage6_schema_exports_required_columns() -> None:
    assert "pose_id" in schemas.IDENTIFIER_COLUMNS
    assert "rank_relevance_label" in schemas.LABEL_COLUMNS
    assert "allowed_for_training" in schemas.LEAKAGE_AUDIT_COLUMNS
    assert "group_usable_for_ranking" in schemas.RANKING_GROUP_COLUMNS
