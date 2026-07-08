from __future__ import annotations

import json

import pandas as pd


def compact_failure_summary(validation: pd.DataFrame, acceptance: pd.DataFrame | None = None) -> str:
    failures = validation["rejection_reason"].fillna("").value_counts().head(5).to_dict() if not validation.empty else {}
    accepted = int(acceptance["accepted_flag"].sum()) if acceptance is not None and not acceptance.empty else 0
    return json.dumps({"dominant_validation_failures": failures, "accepted": accepted}, sort_keys=True)
