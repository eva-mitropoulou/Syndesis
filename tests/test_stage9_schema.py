from egfr_dockingforge.stage9.schemas import TABLE_SCHEMAS


def test_stage9_required_table_schemas_exist():
    required = {
        "analog_seed_scaffolds",
        "edit_sites",
        "transformation_library",
        "analog_candidates",
        "analog_validation",
        "analog_screening_results",
        "analog_acceptance",
        "analog_strategy_benchmark",
    }
    assert required.issubset(TABLE_SCHEMAS)
    for columns in TABLE_SCHEMAS.values():
        assert columns
        assert len(columns) == len(set(columns))
