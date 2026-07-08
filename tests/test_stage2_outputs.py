from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest


PROJECT_ROOT = Path(__file__).resolve().parents[1]


def test_stage2_reference_controls_when_outputs_exist() -> None:
    ensemble_path = PROJECT_ROOT / "data/processed/stage2/receptor_ensemble_v1.parquet"
    exclusion_path = PROJECT_ROOT / "data/processed/stage2/receptor_exclusion_table.parquet"
    if not ensemble_path.exists():
        pytest.skip("Stage 2 ensemble has not been generated.")
    ensemble = pd.read_parquet(ensemble_path)
    excluded = pd.read_parquet(exclusion_path)
    assert (ensemble["stage3_validation_status"] == "pending").all()
    assert ensemble["suggested_docking_box_center"].notna().all()
    assert ensemble["suggested_docking_box_size"].notna().all()
    assert "4ZAU" not in set(ensemble["pdb_id"])
    controls = set(ensemble["pdb_id"]) | set(excluded["pdb_id"])
    assert {"1M17", "1XKK", "4ZAU"}.issubset(controls)

