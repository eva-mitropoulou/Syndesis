from __future__ import annotations

from pathlib import Path
from typing import Any

import pandas as pd

from syndesis.common.io import write_json
from syndesis.stage6.feature_builder import build_pose_features
from syndesis.stage6.interpretability import write_interpretability
from syndesis.stage6.label_builder import build_pose_labels
from syndesis.stage6.leakage_audit import audit_pose_features
from syndesis.stage6.load_stage_inputs import load_stage6_config, load_stage6_inputs, stage6_paths
from syndesis.stage6.model_selection import train_and_evaluate_models
from syndesis.stage6.report_stage6 import write_stage6_report
from syndesis.stage6.splitter import build_ranking_groups, build_splits


def _load(config_path: str | Path) -> tuple[dict[str, Any], dict[str, Path], dict[str, pd.DataFrame]]:
    config = load_stage6_config(config_path)
    paths = stage6_paths(config)
    inputs = load_stage6_inputs(config)
    return config, paths, inputs


def build_pose_features_cli(config_path: str | Path) -> dict[str, Any]:
    config, paths, inputs = _load(config_path)
    features = build_pose_features(inputs, config, paths)
    return {"status": "complete", "rows": int(len(features)), "columns": int(len(features.columns))}


def audit_pose_features_cli(config_path: str | Path) -> dict[str, Any]:
    config, paths, inputs = _load(config_path)
    feature_path = paths["processed"] / "pose_model_features.parquet"
    features = pd.read_parquet(feature_path) if feature_path.exists() else build_pose_features(inputs, config, paths)
    audit = audit_pose_features(features, config, paths)
    return {"status": "complete", "features": int(len(audit)), "trainable": int(audit["allowed_for_training"].sum())}


def build_pose_labels_cli(config_path: str | Path) -> dict[str, Any]:
    config, paths, inputs = _load(config_path)
    labels = build_pose_labels(inputs, config, paths)
    return {"status": "complete", "labels": int(len(labels))}


def split_pose_data_cli(config_path: str | Path) -> dict[str, Any]:
    config, paths, inputs = _load(config_path)
    features = pd.read_parquet(paths["processed"] / "pose_model_features.parquet")
    labels = pd.read_parquet(paths["processed"] / "pose_model_labels.parquet")
    groups = build_ranking_groups(features, labels, paths)
    splits = build_splits(features, inputs, config, paths)
    return {"status": "complete", "groups": int(len(groups)), "split_rows": int(len(splits))}


def train_pose_rankers_cli(config_path: str | Path) -> dict[str, Any]:
    config, paths, _inputs = _load(config_path)
    result = train_and_evaluate_models(
        pd.read_parquet(paths["processed"] / "pose_model_features.parquet"),
        pd.read_parquet(paths["processed"] / "pose_model_labels.parquet"),
        pd.read_parquet(paths["processed"] / "ranking_groups.parquet"),
        pd.read_parquet(paths["processed"] / "model_splits.parquet"),
        pd.read_parquet(paths["processed"] / "feature_leakage_audit.parquet"),
        config,
        paths,
    )
    return {"status": "complete", "artifact": str(result["artifact_path"])}


def train_pose_confidence_classifier_cli(config_path: str | Path) -> dict[str, Any]:
    return train_pose_rankers_cli(config_path)


def calibrate_pose_confidence_cli(config_path: str | Path) -> dict[str, Any]:
    config, paths, _inputs = _load(config_path)
    metrics = pd.read_parquet(paths["processed"] / "calibration_metrics.parquet")
    return {"status": "complete", "calibration_rows": int(len(metrics))}


def evaluate_pose_models_cli(config_path: str | Path) -> dict[str, Any]:
    config, paths, _inputs = _load(config_path)
    metrics = pd.read_parquet(paths["processed"] / "model_metrics.parquet")
    return {"status": "complete", "metrics": int(len(metrics))}


def explain_pose_model_cli(config_path: str | Path) -> dict[str, Any]:
    config, paths, _inputs = _load(config_path)
    fi, ablation = write_interpretability(
        pd.read_parquet(paths["processed"] / "pose_model_features.parquet"),
        pd.read_parquet(paths["processed"] / "pose_model_labels.parquet"),
        pd.read_parquet(paths["processed"] / "feature_leakage_audit.parquet"),
        paths,
    )
    return {"status": "complete", "feature_importance_rows": int(len(fi)), "ablation_rows": int(len(ablation))}


def select_pose_model_cli(config_path: str | Path) -> dict[str, Any]:
    config, paths, _inputs = _load(config_path)
    selection = pd.read_parquet(paths["processed"] / "model_selection_summary.parquet")
    return {"status": "complete", "selected_ranker": str(selection["selected_ranker_model_id"].iloc[0])}


def report_stage6_cli(config_path: str | Path) -> dict[str, Any]:
    config, paths, _inputs = _load(config_path)
    report = write_stage6_report(
        pd.read_parquet(paths["processed"] / "pose_model_features.parquet"),
        pd.read_parquet(paths["processed"] / "pose_model_labels.parquet"),
        pd.read_parquet(paths["processed"] / "ranking_groups.parquet"),
        pd.read_parquet(paths["processed"] / "feature_leakage_audit.parquet"),
        pd.read_parquet(paths["processed"] / "model_splits.parquet"),
        pd.read_parquet(paths["processed"] / "model_metrics.parquet"),
        pd.read_parquet(paths["processed"] / "model_selection_summary.parquet"),
        pd.read_parquet(paths["processed"] / "feature_importance.parquet"),
        pd.read_parquet(paths["processed"] / "feature_group_ablation.parquet"),
        paths["reports"] / "06_pose_reranking_confidence.html",
    )
    return {"status": "complete", "report": str(report)}


def run_stage6_all(config_path: str | Path) -> dict[str, Any]:
    config, paths, inputs = _load(config_path)
    features = build_pose_features(inputs, config, paths)
    audit = audit_pose_features(features, config, paths)
    labels = build_pose_labels(inputs, config, paths)
    groups = build_ranking_groups(features, labels, paths)
    splits = build_splits(features, inputs, config, paths)
    result = train_and_evaluate_models(features, labels, groups, splits, audit, config, paths)
    fi, ablation = write_interpretability(features, labels, audit, paths)
    report = write_stage6_report(features, labels, groups, audit, splits, result["metrics"], result["selection"], fi, ablation, paths["reports"] / "06_pose_reranking_confidence.html")
    payload = {
        "status": "complete",
        "features": int(len(features)),
        "labels": int(len(labels)),
        "groups": int(len(groups)),
        "report": str(report),
        "model_artifact": str(result["artifact_path"]),
    }
    write_json(paths["processed"] / "stage6_summary.json", payload)
    return payload
