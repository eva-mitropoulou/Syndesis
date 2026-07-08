from __future__ import annotations

from pathlib import Path

import pandas as pd

from egfr_dockingforge.stage1.schemas import BENCHMARK_COLUMNS, QUALITY_TIERS, empty_benchmark_frame


PROJECT_ROOT = Path(__file__).resolve().parents[1]


def test_benchmark_schema_has_required_columns() -> None:
    frame = empty_benchmark_frame()
    required = {
        "complex_id",
        "pdb_id",
        "ligand_comp_id",
        "native_ligand_sdf_path",
        "receptor_clean_path",
        "quality_tier",
        "quality_score",
        "include_in_stage1_benchmark",
        "exclusion_reason",
        "warnings_json",
    }
    assert required.issubset(set(frame.columns))
    assert set(BENCHMARK_COLUMNS).issubset(set(frame.columns))
    assert isinstance(frame, pd.DataFrame)


def test_rejected_table_target_path_is_defined() -> None:
    config_path = PROJECT_ROOT / "configs/stage1_cocrystal_benchmark.yaml"
    text = config_path.read_text(encoding="utf-8")
    assert "data/interim/stage1" in text
    assert "data/processed/stage1" in text


def test_quality_tier_values_are_defined() -> None:
    assert QUALITY_TIERS == {"Tier A", "Tier B", "Tier C", "Rejected"}

