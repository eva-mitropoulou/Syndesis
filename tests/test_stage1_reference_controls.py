from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest


PROJECT_ROOT = Path(__file__).resolve().parents[1]


def test_stage1_reference_controls_when_benchmark_exists() -> None:
    benchmark_path = PROJECT_ROOT / "data/processed/stage1/egfr_cocrystal_benchmark.csv"
    if not benchmark_path.exists():
        pytest.skip("Stage 1 benchmark has not been generated in this checkout.")

    benchmark = pd.read_csv(benchmark_path)
    controls = benchmark[benchmark["pdb_id"].isin(["1M17", "1XKK", "4ZAU"])]
    assert {"1M17", "1XKK", "4ZAU"}.issubset(set(controls["pdb_id"]))

    positives = controls[controls["pdb_id"].isin(["1M17", "1XKK"])]
    assert not positives.empty
    assert positives["include_in_stage1_benchmark"].isin([True]).all()

    excluded = controls[controls["pdb_id"] == "4ZAU"]
    assert not excluded.empty
    assert excluded["include_in_stage1_benchmark"].isin([False]).all()
    assert excluded["exclusion_reason"].fillna("").str.contains("Covalent", case=False).any()

