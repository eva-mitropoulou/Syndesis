from syndesis.stage10.metric_definitions import accepted_pre_md


def test_hard_scope_failure_cannot_be_accepted():
    row = {
        "valid_molecule_flag": True,
        "unique_flag": True,
        "hard_scope_pass": False,
        "covalent_warhead_flag": False,
        "reactive_flag": False,
        "medchem_risk_score": 0.0,
        "binding_mode_preserved_flag": True,
        "best_pose_confidence": 0.9,
        "delta_candidate_score": 0.1,
        "delta_ligand_efficiency": 0.1,
    }
    assert not accepted_pre_md(row)
