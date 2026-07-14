from __future__ import annotations

from pathlib import Path
from typing import Any

import pandas as pd

from syndesis.common.io import write_json
from syndesis.stage8.candidate_aggregation import aggregate_candidates, pose_decisions
from syndesis.stage8.failure_analysis import summarize_failures
from syndesis.stage8.load_stage_inputs import load_stage8_config, load_stage8_inputs, stage8_paths
from syndesis.stage8.pose_confidence_apply import apply_pose_confidence
from syndesis.stage8.report_stage8 import write_stage8_report
from syndesis.stage8.score_normalization import normalize_scores, primary_triage
from syndesis.stage8.screening_docking import run_candidate_docking
from syndesis.stage8.screening_interactions import compute_screening_interactions, pose_sanity
from syndesis.stage8.screening_rescoring import rescore_screening_poses
from syndesis.stage8.screening_task_matrix import build_screening_task_matrix, prepare_manifest
from syndesis.stage8.source_aware_selection import rank_candidates


def _load(config_path: str | Path) -> tuple[dict[str, Any], dict[str, Path], dict[str, pd.DataFrame]]:
    config = load_stage8_config(config_path)
    paths = stage8_paths(config)
    inputs = load_stage8_inputs(config)
    return config, paths, inputs


def build_screening_task_matrix_cli(config_path: str | Path) -> dict:
    config, paths, inputs = _load(config_path)
    manifest = prepare_manifest(inputs["stage7_screening_input"], inputs["candidate_library_master"], config, paths)
    tasks = build_screening_task_matrix(manifest, inputs["receptor_ensemble"], config, paths)
    return {"status": "complete", "manifest_rows": len(manifest), "tasks": len(tasks)}


def run_candidate_docking_cli(config_path: str | Path) -> dict:
    config, paths, _inputs = _load(config_path)
    tasks = pd.read_parquet(paths["processed"] / "screening_task_matrix.parquet")
    runs, poses = run_candidate_docking(tasks, config, paths)
    triage = primary_triage(poses, config, paths)
    return {"status": "complete", "runs": len(runs), "poses": len(poses), "triaged": len(triage)}


def rescore_screening_poses_cli(config_path: str | Path) -> dict:
    config, paths, _inputs = _load(config_path)
    scores = rescore_screening_poses(pd.read_parquet(paths["processed"] / "screening_docked_poses.parquet"), pd.read_parquet(paths["processed"] / "primary_docking_triage.parquet"), pd.read_parquet(paths["processed"] / "screening_task_matrix.parquet"), config, paths)
    return {"status": "complete", "scores": len(scores)}


def compute_screening_interactions_cli(config_path: str | Path) -> dict:
    config, paths, inputs = _load(config_path)
    poses = pd.read_parquet(paths["processed"] / "screening_docked_poses.parquet")
    sanity = pose_sanity(poses, paths)
    long, features = compute_screening_interactions(poses, pd.read_parquet(paths["processed"] / "screening_gnina_scores.parquet"), pd.read_parquet(paths["processed"] / "screening_task_matrix.parquet"), pd.read_parquet(paths["processed"] / "screening_candidate_manifest.parquet"), inputs["residue_map"], inputs["key_interactions"], inputs["native_fingerprints"], config, paths)
    return {"status": "complete", "sanity": len(sanity), "interactions": len(long), "features": len(features)}


def apply_pose_confidence_cli(config_path: str | Path) -> dict:
    config, paths, _inputs = _load(config_path)
    _, conf = apply_pose_confidence(pd.read_parquet(paths["processed"] / "screening_docked_poses.parquet"), pd.read_parquet(paths["processed"] / "screening_gnina_scores.parquet"), pd.read_parquet(paths["processed"] / "screening_pose_sanity.parquet"), pd.read_parquet(paths["processed"] / "screening_interaction_features.parquet"), pd.read_parquet(paths["processed"] / "screening_task_matrix.parquet"), pd.read_parquet(paths["processed"] / "screening_candidate_manifest.parquet"), config, paths)
    normalize_scores(conf, paths)
    return {"status": "complete", "confidence_rows": len(conf)}


def aggregate_candidate_scores_cli(config_path: str | Path) -> dict:
    config, paths, inputs = _load(config_path)
    conf = pd.read_parquet(paths["processed"] / "screening_pose_confidence.parquet")
    norm = pd.read_parquet(paths["processed"] / "normalized_screening_scores.parquet")
    manifest = pd.read_parquet(paths["processed"] / "screening_candidate_manifest.parquet")
    decisions = pose_decisions(conf, norm, manifest, config, paths)
    agg = aggregate_candidates(conf, decisions, manifest, inputs["candidate_library_master"], paths)
    return {"status": "complete", "aggregate_rows": len(agg)}


def select_ranked_candidates_cli(config_path: str | Path) -> dict:
    config, paths, _inputs = _load(config_path)
    ranked, diag = rank_candidates(pd.read_parquet(paths["processed"] / "candidate_aggregate_scores.parquet"), paths)
    summarize_failures(paths)
    return {"status": "complete", "ranked": len(ranked), "diagnostics": len(diag)}


def report_stage8_cli(config_path: str | Path) -> dict:
    config, paths, _inputs = _load(config_path)
    report = write_stage8_report(paths)
    return {"status": "complete", "report": str(report)}


def run_stage8_all(config_path: str | Path) -> dict:
    build_screening_task_matrix_cli(config_path)
    run_candidate_docking_cli(config_path)
    rescore_screening_poses_cli(config_path)
    compute_screening_interactions_cli(config_path)
    apply_pose_confidence_cli(config_path)
    aggregate_candidate_scores_cli(config_path)
    select_ranked_candidates_cli(config_path)
    summary = report_stage8_cli(config_path)
    config, paths, _inputs = _load(config_path)
    write_json(paths["processed"] / "stage8_summary.json", summary)
    return summary
