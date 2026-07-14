from __future__ import annotations

import pandas as pd

from syndesis.stage8.candidate_aggregation import pose_decisions


def test_pose_decision_does_not_use_docking_score_alone(tmp_path) -> None:
    conf = pd.DataFrame({"screening_pose_id": ["p1"], "molecule_id": ["m1"], "target_receptor_id": ["r"], "receptor_state": ["active"], "sanity_status": ["pass"], "pose_confidence_probability": [0.1], "key_interaction_recall_consensus": [0.9], "warnings_json": ["[]"]})
    norm = pd.DataFrame({"screening_pose_id": ["p1"], "molecule_id": ["m1"], "target_receptor_id": ["r"], "receptor_state": ["active"]})
    out = pose_decisions(conf, norm, pd.DataFrame(), {"triage": {"confidence_threshold": 0.3, "interaction_f1_threshold": 0.2}}, {"processed": tmp_path})
    assert out["pose_decision_label"].iloc[0] == "low_confidence"
