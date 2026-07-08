from __future__ import annotations

import pandas as pd

from egfr_dockingforge.stage8.score_normalization import percentile


def test_percentile_direction() -> None:
    s = pd.Series([1.0, 2.0, 3.0])
    assert percentile(s, higher_is_better=True).iloc[-1] == 1.0
    assert percentile(s, higher_is_better=False).iloc[0] == 1.0
