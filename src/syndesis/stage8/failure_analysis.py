from __future__ import annotations

import pandas as pd


def summarize_failures(paths: dict) -> pd.DataFrame:
    runs = pd.read_parquet(paths["processed"] / "screening_docking_runs.parquet")
    rows = runs.groupby("status").size().reset_index(name="count") if not runs.empty else pd.DataFrame(columns=["status", "count"])
    rows.to_parquet(paths["processed"] / "screening_failure_summary.parquet", index=False)
    rows.to_csv(paths["processed"] / "screening_failure_summary.csv", index=False)
    return rows
