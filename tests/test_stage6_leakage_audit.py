from __future__ import annotations

import pytest

from egfr_dockingforge.stage6.leakage_audit import assert_no_leakage, classify_feature


def test_forbidden_native_features_are_rejected() -> None:
    with pytest.raises(RuntimeError):
        assert_no_leakage(["cnnscore", "rmsd_symmetry_corrected"])


def test_classify_feature_marks_native_key_recall_as_drop() -> None:
    group, train, deploy, risk, _reason, action = classify_feature(
        "key_interaction_recall_native", {"features": {"forbidden_patterns": ["native"]}}
    )
    assert group == "label_or_native_diagnostic"
    assert train is False
    assert deploy is False
    assert risk == "high"
    assert action == "drop_for_training"
