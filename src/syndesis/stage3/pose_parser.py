from __future__ import annotations

import pandas as pd

from syndesis.stage3.schemas import DOCKED_POSE_COLUMNS


def empty_pose_table() -> pd.DataFrame:
    return pd.DataFrame(columns=DOCKED_POSE_COLUMNS)

