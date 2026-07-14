from __future__ import annotations

import pandas as pd


def deduplicate_analogs(candidates: pd.DataFrame) -> pd.DataFrame:
    out = candidates.copy()
    dupes = out.duplicated(["standard_smiles"], keep="first")
    out.loc[dupes, "uniqueness_status"] = "duplicate_global"
    return out.drop_duplicates(["standard_smiles"], keep="first").reset_index(drop=True)
