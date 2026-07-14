from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from egfr_dockingforge.common.io import write_table


def bootstrap_ci(values: list[float], iterations: int = 1000, seed: int = 1010) -> tuple[float, float]:
    if not values:
        return 0.0, 0.0
    rng = np.random.default_rng(seed)
    arr = np.array(values, dtype=float)
    boots = [float(rng.choice(arr, size=len(arr), replace=True).mean()) for _ in range(iterations)]
    return float(np.percentile(boots, 2.5)), float(np.percentile(boots, 97.5))


def benjamini_hochberg(p_values: list[float]) -> list[float]:
    n = len(p_values)
    if n == 0:
        return []
    order = np.argsort(p_values)
    corrected = np.empty(n)
    prev = 1.0
    for rank, idx in reversed(list(enumerate(order, start=1))):
        val = min(prev, p_values[idx] * n / rank)
        corrected[idx] = val
        prev = val
    return corrected.tolist()


def paired_permutation_pvalue(a: list[float], b: list[float]) -> float:
    if len(a) != len(b) or not a:
        return 1.0
    diffs = np.array(a) - np.array(b)
    observed = abs(diffs.mean())
    if len(diffs) > 12:
        rng = np.random.default_rng(1010)
        samples = [abs((diffs * rng.choice([-1, 1], size=len(diffs))).mean()) for _ in range(2000)]
    else:
        samples = []
        for mask in range(2 ** len(diffs)):
            signs = np.array([1 if (mask >> i) & 1 else -1 for i in range(len(diffs))])
            samples.append(abs((diffs * signs).mean()))
    return float((np.array(samples) >= observed).mean()) if samples else 1.0


def run_statistical_comparisons(seed_metrics: pd.DataFrame, config: dict[str, Any], paths: dict[str, Path]) -> pd.DataFrame:
    primary = config["benchmark"]["primary_method"]
    methods = [
        "random_analog_enumeration","docking_score_only_optimization","gnina_only_optimization","rdkit_rule_based",
    ]
    rows = []
    for method in methods:
        a = seed_metrics[seed_metrics["strategy_name"].eq(primary)].set_index("seed_id")
        b = seed_metrics[seed_metrics["strategy_name"].eq(method)].set_index("seed_id")
        seeds = sorted(set(a.index) & set(b.index))
        avals = [float(a.loc[s, "accepted_rate_pre_md"]) for s in seeds]
        bvals = [float(b.loc[s, "accepted_rate_pre_md"]) for s in seeds]
        ci = bootstrap_ci([x - y for x, y in zip(avals, bvals)], int(config["benchmark"]["bootstrap_iterations"]), int(config["benchmark"]["random_seed"]))
        p = paired_permutation_pvalue(avals, bvals)
        rows.append({
            "comparison_id": f"{primary}_vs_{method}",
            "metric_name": "accepted_rate_pre_md",
            "method_a": primary,
            "method_b": method,
            "mean_a": float(np.mean(avals)) if avals else 0.0,
            "mean_b": float(np.mean(bvals)) if bvals else 0.0,
            "median_a": float(np.median(avals)) if avals else 0.0,
            "median_b": float(np.median(bvals)) if bvals else 0.0,
            "delta_mean": float(np.mean(avals) - np.mean(bvals)) if avals and bvals else 0.0,
            "effect_size": float(np.mean(avals) - np.mean(bvals)) if avals and bvals else 0.0,
            "ci_low": ci[0],
            "ci_high": ci[1],
            "p_value": p,
            "p_value_corrected": p,
            "test_name": "paired_permutation_over_seed_scaffolds",
            "num_seed_pairs": len(seeds),
            "interpretation": "exploratory_small_seed_count" if len(seeds) < 5 else "paired_budget_controlled_comparison",
            "warnings_json": json.dumps(["small_seed_count"] if len(seeds) < 5 else []),
        })
    corrected = benjamini_hochberg([row["p_value"] for row in rows])
    for row, corr in zip(rows, corrected):
        row["p_value_corrected"] = corr
    out = pd.DataFrame(rows)
    write_table(paths["processed"] / "statistical_comparisons.parquet", out)
    write_table(paths["processed"] / "statistical_comparisons.csv", out)
    return out
