from __future__ import annotations

from pathlib import Path
from typing import Any

import pandas as pd

from egfr_dockingforge.common.io import write_json
from egfr_dockingforge.stage10.ablation_tables import build_ablation_summary
from egfr_dockingforge.stage10.acceptance_metrics import build_master_table, compute_seed_and_strategy_metrics
from egfr_dockingforge.stage10.benchmark_manifest import build_ablation_manifest, build_budget_audit
from egfr_dockingforge.stage10.compute_cost_metrics import compute_cost_metrics
from egfr_dockingforge.stage10.diversity_metrics import compute_diversity_novelty
from egfr_dockingforge.stage10.load_stage_inputs import load_stage10_config, load_stage10_inputs, stage10_paths
from egfr_dockingforge.stage10.plots import make_ablation_plots
from egfr_dockingforge.stage10.report_stage10 import write_stage10_report
from egfr_dockingforge.stage10.score_hacking_metrics import compute_score_hacking_cases
from egfr_dockingforge.stage10.statistics import run_statistical_comparisons


def _load(config_path: str | Path) -> tuple[dict[str, Any], dict[str, Path], dict[str, pd.DataFrame | None]]:
    config = load_stage10_config(config_path)
    paths = stage10_paths(config)
    inputs = load_stage10_inputs(config)
    return config, paths, inputs


def build_ablation_manifest_cli(config_path: str | Path) -> dict:
    config, paths, inputs = _load(config_path)
    manifest = build_ablation_manifest(inputs, config, paths)
    master = build_master_table(inputs, manifest, config, paths)
    budget = build_budget_audit(master, manifest, inputs, config, paths)
    return {"status": "complete", "strategies": len(manifest), "budget_rows": len(budget)}


def compute_analog_benchmark_metrics_cli(config_path: str | Path) -> dict:
    config, paths, inputs = _load(config_path)
    manifest = pd.read_parquet(paths["processed"] / "ablation_strategy_manifest.parquet")
    master = pd.read_parquet(paths["processed"] / "analog_benchmark_master.parquet")
    seed_metrics, strategy_metrics = compute_seed_and_strategy_metrics(master, manifest, paths)
    compute_diversity_novelty(master, paths)
    budget = pd.read_parquet(paths["processed"] / "strategy_budget_audit.parquet")
    compute_cost_metrics(strategy_metrics, budget, paths)
    return {"status": "complete", "seed_rows": len(seed_metrics), "strategy_rows": len(strategy_metrics)}


def compute_score_hacking_metrics_cli(config_path: str | Path) -> dict:
    config, paths, inputs = _load(config_path)
    master = pd.read_parquet(paths["processed"] / "analog_benchmark_master.parquet")
    cases = compute_score_hacking_cases(master, paths)
    return {"status": "complete", "cases": len(cases)}


def run_ablation_statistics_cli(config_path: str | Path) -> dict:
    config, paths, inputs = _load(config_path)
    seed_metrics = pd.read_parquet(paths["processed"] / "seed_strategy_metrics.parquet")
    strategy_metrics = pd.read_parquet(paths["processed"] / "strategy_metrics.parquet")
    stats = run_statistical_comparisons(seed_metrics, config, paths)
    ablations = build_ablation_summary(strategy_metrics, paths, inputs)
    return {"status": "complete", "comparisons": len(stats), "ablations": len(ablations)}


def make_ablation_plots_cli(config_path: str | Path) -> dict:
    config, paths, inputs = _load(config_path)
    strategy_metrics = pd.read_parquet(paths["processed"] / "strategy_metrics.parquet")
    figs = make_ablation_plots(strategy_metrics, paths)
    return {"status": "complete", "figures": len(figs)}


def report_stage10_cli(config_path: str | Path) -> dict:
    config, paths, inputs = _load(config_path)
    report = write_stage10_report(paths)
    return {"status": "complete", "report": str(report)}


def run_stage10_all(config_path: str | Path) -> dict:
    build_ablation_manifest_cli(config_path)
    compute_analog_benchmark_metrics_cli(config_path)
    compute_score_hacking_metrics_cli(config_path)
    run_ablation_statistics_cli(config_path)
    make_ablation_plots_cli(config_path)
    summary = report_stage10_cli(config_path)
    config, paths, _ = _load(config_path)
    write_json(paths["processed"] / "stage10_summary.json", summary)
    return summary
