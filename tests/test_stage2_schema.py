from __future__ import annotations

from egfr_dockingforge.stage2.schemas import ENSEMBLE_COLUMNS, FEATURE_COLUMNS, empty_ensemble_frame, empty_features_frame


def test_stage2_feature_schema_contains_required_columns() -> None:
    frame = empty_features_frame()
    required = {
        "receptor_id",
        "complex_id",
        "pdb_id",
        "state_stratum",
        "lys745_glu762_nz_oe_min_distance",
        "native_ligand_centroid_x",
        "receptor_preselection_score",
        "warnings_json",
    }
    assert required.issubset(frame.columns)
    assert set(FEATURE_COLUMNS).issubset(frame.columns)


def test_stage2_ensemble_schema_contains_required_columns() -> None:
    frame = empty_ensemble_frame()
    required = {
        "receptor_id",
        "selected_flag",
        "selected_role",
        "suggested_docking_box_center",
        "suggested_docking_box_size",
        "stage3_validation_status",
    }
    assert required.issubset(frame.columns)
    assert set(ENSEMBLE_COLUMNS).issubset(frame.columns)

