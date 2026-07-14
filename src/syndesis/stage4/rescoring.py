from __future__ import annotations

from pathlib import Path
from typing import Any

import pandas as pd

from syndesis.common.io import write_json
from syndesis.stage4.gnina_runner import run_gnina_rescoring
from syndesis.stage4.load_stage3 import load_stage3_inputs, load_stage4_config, stage4_paths
from syndesis.stage4.report_stage4 import write_stage4_report
from syndesis.stage4.score_diagnostics import (
    build_pose_score_table,
    failure_table,
    rescoring_task_metrics,
    scorer_summary,
)
from syndesis.stage4.scoring_task_matrix import build_rescoring_tasks
from syndesis.stage4.vina_rescore_runner import build_empirical_scores


def build_rescoring_task_matrix(config_path: str | Path) -> dict[str, Any]:
    config = load_stage4_config(config_path)
    paths = stage4_paths(config)
    inputs = load_stage3_inputs(config)
    tasks = build_rescoring_tasks(inputs, config, paths)
    return {"status": "complete", "tasks": int(len(tasks)), "ready": int(tasks["task_status"].eq("ready").sum())}


def run_gnina_stage4(config_path: str | Path) -> dict[str, Any]:
    config = load_stage4_config(config_path)
    paths = stage4_paths(config)
    task_path = paths["processed"] / "rescoring_task_matrix.parquet"
    if not task_path.exists():
        build_rescoring_task_matrix(config_path)
    tasks = pd.read_parquet(task_path)
    raw, scores = run_gnina_rescoring(tasks, config, paths)
    return {"status": "complete", "runs": int(len(raw)), "scores": int(len(scores))}


def run_empirical_rescoring(config_path: str | Path) -> dict[str, Any]:
    config = load_stage4_config(config_path)
    paths = stage4_paths(config)
    tasks = pd.read_parquet(paths["processed"] / "rescoring_task_matrix.parquet")
    gnina_path = paths["processed"] / "gnina_scores.parquet"
    gnina = pd.read_parquet(gnina_path) if gnina_path.exists() else pd.DataFrame()
    empirical = build_empirical_scores(tasks, gnina, paths)
    return {"status": "complete", "scores": int(len(empirical))}


def parse_rescoring_outputs(config_path: str | Path) -> dict[str, Any]:
    config = load_stage4_config(config_path)
    paths = stage4_paths(config)
    tasks = pd.read_parquet(paths["processed"] / "rescoring_task_matrix.parquet")
    labels = load_stage3_inputs(config)["labels"]
    gnina = pd.read_parquet(paths["processed"] / "gnina_scores.parquet")
    empirical_path = paths["processed"] / "empirical_scores.parquet"
    empirical = pd.read_parquet(empirical_path) if empirical_path.exists() else build_empirical_scores(tasks, gnina, paths)
    table = build_pose_score_table(tasks, labels, empirical, gnina, config, paths)
    return {"status": "complete", "pose_scores": int(len(table))}


def diagnose_rescoring(config_path: str | Path) -> dict[str, Any]:
    config = load_stage4_config(config_path)
    paths = stage4_paths(config)
    score_table = pd.read_parquet(paths["processed"] / "pose_score_table.parquet")
    raw_path = paths["processed"] / "gnina_raw_runs.parquet"
    raw = pd.read_parquet(raw_path) if raw_path.exists() else pd.DataFrame()
    metrics = rescoring_task_metrics(score_table, config, paths)
    summary = scorer_summary(score_table, metrics, raw, config, paths)
    failures = failure_table(score_table, paths)
    return {"status": "complete", "task_metrics": int(len(metrics)), "summary_rows": int(len(summary)), "failures": int(len(failures))}


def report_stage4(config_path: str | Path) -> dict[str, Any]:
    config = load_stage4_config(config_path)
    paths = stage4_paths(config)
    tasks = pd.read_parquet(paths["processed"] / "rescoring_task_matrix.parquet")
    gnina = pd.read_parquet(paths["processed"] / "gnina_scores.parquet")
    metrics = pd.read_parquet(paths["processed"] / "rescoring_task_metrics.parquet")
    summary = pd.read_parquet(paths["processed"] / "scorer_comparison_summary.parquet")
    failures = pd.read_parquet(paths["processed"] / "rescoring_failures.parquet")
    raw_path = paths["processed"] / "gnina_raw_runs.parquet"
    raw = pd.read_parquet(raw_path) if raw_path.exists() else pd.DataFrame()
    report = write_stage4_report(tasks, gnina, metrics, summary, failures, raw, paths["reports"] / "04_ml_rescoring.html")
    payload = {"status": "complete", "tasks": int(len(tasks)), "gnina_scores": int(len(gnina)), "report": str(report)}
    write_json(paths["processed"] / "stage4_summary.json", payload)
    return payload


def build_stage4_all(config_path: str | Path) -> dict[str, Any]:
    build_rescoring_task_matrix(config_path)
    run_gnina_stage4(config_path)
    run_empirical_rescoring(config_path)
    parse_rescoring_outputs(config_path)
    diagnose_rescoring(config_path)
    return report_stage4(config_path)

