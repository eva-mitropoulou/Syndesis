from __future__ import annotations

import pandas as pd

from egfr_dockingforge.stage8.screening_task_matrix import _vec


def test_vec_parses_box_lists() -> None:
    assert _vec("[1.0, 2.0, 3.0]") == [1.0, 2.0, 3.0]
