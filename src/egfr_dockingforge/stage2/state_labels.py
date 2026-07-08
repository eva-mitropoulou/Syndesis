from __future__ import annotations

from typing import Any

import pandas as pd


def normalize_state_label(row: pd.Series) -> dict[str, Any]:
    kincore = row.get("kincore_state")
    dfg = row.get("dfg_state")
    activation = row.get("activation_loop_state")
    if isinstance(dfg, str) and dfg.lower() in {"dfgout", "dfg-out", "out"}:
        stratum = "dfgout_or_typeII_like"
        activity = "inactive-like"
    elif isinstance(kincore, str) and "active" in kincore.lower() and "inactive" not in kincore.lower():
        stratum = "active_like"
        activity = "active-like"
    elif isinstance(kincore, str) and "inactive" in kincore.lower():
        stratum = "inactive_like"
        activity = "inactive-like"
    elif isinstance(activation, str) and "inactive" in activation.lower():
        stratum = "inactive_like"
        activity = "inactive-like"
    else:
        stratum = "unknown_state"
        activity = None
    return {
        "kincore_activity_label": activity,
        "state_stratum": stratum,
        "dfg_dihedral_cluster": None,
    }
