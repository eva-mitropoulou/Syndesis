from egfr_dockingforge.stage10.metric_definitions import score_hacking


def test_score_improves_and_pose_confidence_worsens_is_hacking():
    assert score_hacking(
        {
            "delta_gnina_cnnscore": 0.2,
            "delta_pose_confidence": -0.1,
            "delta_key_interaction_recall": 0.0,
            "binding_mode_preserved_flag": True,
            "medchem_risk_score": 0.0,
        }
    )


def test_score_improves_and_binding_mode_breaks_is_hacking():
    assert score_hacking({"delta_gnina_cnnscore": 0.2, "binding_mode_preserved_flag": False, "medchem_risk_score": 0.0})
