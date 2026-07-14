from syndesis.stage11.schemas import TABLE_SCHEMAS


def test_stage11_schemas_exist():
    required = {"md_candidate_manifest","ligand_parameterization","ligand_parameterization_report","system_builds","md_runs","equilibration_qc","md_metrics","md_interaction_persistence","md_binding_mode_persistence","md_pose_stability_labels","md_candidate_summary","stage10_post_md_input"}
    assert required.issubset(TABLE_SCHEMAS)
    assert all(len(cols) == len(set(cols)) for cols in TABLE_SCHEMAS.values())
