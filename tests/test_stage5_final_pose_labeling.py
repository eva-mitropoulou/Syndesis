from __future__ import annotations

import pandas as pd

from syndesis.stage5.final_pose_labeling import classify_final_pose


def test_stage5_high_confidence_requires_rmsd_sanity_and_recovery() -> None:
    row = pd.Series({"rmsd_symmetry_corrected": 1.5, "sanity_status": "pass", "interaction_recovery_label": "high_recovery", "key_interaction_recall_native": 1.0})
    label, _confidence, _reason = classify_final_pose(row, {})
    assert label == "high_confidence_native_like"


def test_stage5_rmsd_good_interactions_poor_is_separate() -> None:
    row = pd.Series({"rmsd_symmetry_corrected": 1.5, "sanity_status": "pass", "interaction_recovery_label": "poor_recovery", "key_interaction_recall_native": 0.0})
    label, _confidence, _reason = classify_final_pose(row, {})
    assert label == "rmsd_good_interactions_poor"


def test_stage5_invalid_pose_cannot_be_high_confidence() -> None:
    row = pd.Series({"rmsd_symmetry_corrected": 1.0, "sanity_status": "severe_fail", "interaction_recovery_label": "high_recovery"})
    label, _confidence, _reason = classify_final_pose(row, {})
    assert label == "invalid_pose"

