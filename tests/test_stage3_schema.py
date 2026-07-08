from __future__ import annotations

from egfr_dockingforge.stage3.schemas import LABEL_COLUMNS, RECEPTOR_PREP_COLUMNS, TASK_COLUMNS


def test_stage3_required_schema_columns_exist() -> None:
    assert {"receptor_id", "prepared_receptor_file", "preparation_status"}.issubset(RECEPTOR_PREP_COLUMNS)
    assert {"docking_task_id", "task_type", "seed", "reference_pose_transform_matrix_path"}.issubset(TASK_COLUMNS)
    assert {"pose_id", "stage3_pose_label", "interaction_recovery_status"}.issubset(LABEL_COLUMNS)

