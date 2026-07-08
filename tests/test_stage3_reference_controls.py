from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest


PROJECT_ROOT = Path(__file__).resolve().parents[1]


def test_stage3_controls_when_task_matrix_exists() -> None:
    path = PROJECT_ROOT / "data/processed/stage3/docking_task_matrix.parquet"
    if not path.exists():
        pytest.skip("Stage 3 task matrix has not been generated.")
    tasks = pd.read_parquet(path)
    joined = " ".join(tasks["native_receptor_id"].tolist() + tasks["target_receptor_id"].tolist()).lower()
    assert "1m17" in joined
    assert "1xkk" in joined
    assert "4zau" not in joined
    assert tasks["seed"].notna().all()

