from __future__ import annotations

import pandas as pd

from egfr_dockingforge.stage2.state_labels import normalize_state_label


def test_active_like_state_label_normalization() -> None:
    labels = normalize_state_label(pd.Series({"kincore_state": "active-like", "dfg_state": "DFGin"}))
    assert labels["state_stratum"] == "active_like"
    assert labels["kincore_activity_label"] == "active-like"


def test_missing_state_label_becomes_unknown() -> None:
    labels = normalize_state_label(pd.Series({"kincore_state": None, "dfg_state": None}))
    assert labels["state_stratum"] == "unknown_state"

