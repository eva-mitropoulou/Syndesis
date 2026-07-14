from __future__ import annotations

from pathlib import Path
from typing import Any

import pandas as pd

from syndesis.common.io import load_yaml, project_root, resolve_path, write_table


def load_stage2_config(config_path: str | Path) -> dict[str, Any]:
    return load_yaml(resolve_path(config_path, project_root()))


def stage2_paths(config: dict[str, Any]) -> dict[str, Path]:
    root = project_root()
    return {key: resolve_path(value, root) for key, value in config["paths"].items()}


def read_table(path: str | Path) -> pd.DataFrame:
    path = Path(path)
    if path.suffix == ".parquet":
        return pd.read_parquet(path)
    return pd.read_csv(path)


def is_true(value: Any) -> bool:
    return False if pd.isna(value) else bool(value)


def load_stage1_benchmark(config: dict[str, Any]) -> pd.DataFrame:
    path = resolve_path(config["inputs"]["stage1_benchmark"], project_root())
    if not path.exists() and path.with_suffix(".csv").exists():
        path = path.with_suffix(".csv")
    return read_table(path)


def stage2_filter_reason(row: pd.Series, config: dict[str, Any]) -> str:
    filters = config["filters"]
    reasons: list[str] = []
    if is_true(row.get("include_in_stage1_benchmark")) is not True:
        reasons.append("Stage 1 did not include this complex.")
    if filters.get("exclude_covalent", True) and is_true(row.get("covalent_flag")):
        reasons.append("Covalent receptor complex excluded by Stage 2.")
    if filters.get("exclude_allosteric", True) and is_true(row.get("allosteric_flag")):
        reasons.append("Allosteric receptor complex excluded by Stage 2.")
    if is_true(row.get("atp_site_flag")) is not True:
        reasons.append("ATP-site flag is not true.")
    for path_col in ("receptor_clean_path", "native_ligand_sdf_path"):
        value = row.get(path_col)
        if not isinstance(value, str) or not Path(value).exists():
            reasons.append(f"Required file missing: {path_col}.")
    allowed = set(filters.get("allowed_quality_tiers", ["Tier A", "Tier B"]))
    tier = row.get("quality_tier")
    if tier not in allowed:
        if not (tier == "Tier C" and filters.get("allow_tier_c_state_diversity", False)):
            reasons.append(f"Quality tier {tier} is not allowed.")
    resolution = row.get("resolution_angstrom")
    if pd.notna(resolution) and float(resolution) > float(filters.get("max_resolution_angstrom", 3.0)):
        reasons.append("Resolution exceeds Stage 2 maximum.")
    completeness = row.get("active_site_completeness_score")
    if pd.isna(completeness) or float(completeness) < float(filters.get("min_active_site_completeness_score", 0.75)):
        reasons.append("Active-site completeness below Stage 2 threshold.")
    if not filters.get("allow_mutants", True) and is_true(row.get("mutation_flag")):
        reasons.append("Mutant receptors are disabled by config.")
    return "; ".join(reasons)


def filter_stage1_candidates(benchmark: pd.DataFrame, config: dict[str, Any], exclusion_path: Path) -> pd.DataFrame:
    frame = benchmark.copy()
    frame["stage2_exclusion_reason"] = frame.apply(lambda row: stage2_filter_reason(row, config), axis=1)
    frame["stage2_pass_filter"] = frame["stage2_exclusion_reason"].eq("")
    excluded = frame[~frame["stage2_pass_filter"]].copy()
    write_table(exclusion_path.with_suffix(".csv"), excluded)
    write_table(exclusion_path, excluded)
    return frame[frame["stage2_pass_filter"]].copy()
