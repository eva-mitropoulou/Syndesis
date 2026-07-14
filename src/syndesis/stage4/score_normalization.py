from __future__ import annotations

from typing import Any

import pandas as pd


def sort_by_scorer(frame: pd.DataFrame, scorer: str, directions: dict[str, Any]) -> pd.DataFrame:
    direction = str(directions.get(scorer, "lower")).lower()
    ascending = direction != "higher"
    return frame.assign(_score=pd.to_numeric(frame[scorer], errors="coerce")).sort_values(
        ["_score", "original_pose_rank"],
        ascending=[ascending, True],
        na_position="last",
    ).drop(columns=["_score"])

