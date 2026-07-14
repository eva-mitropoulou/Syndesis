from __future__ import annotations

import pandas as pd


def medchem_rejection_rate(master: pd.DataFrame) -> float:
    if master.empty:
        return 0.0
    return float(((master["hard_scope_pass"] == False) | (master["medchem_risk_score"] > 0.4)).mean())
