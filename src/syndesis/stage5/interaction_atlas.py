from __future__ import annotations

from pathlib import Path
from typing import Any

import pandas as pd

from syndesis.common.io import write_json
from syndesis.stage5.binding_mode_clustering import cluster_binding_modes
from syndesis.stage5.complex_builder import build_interaction_complexes
from syndesis.stage5.docked_pose_fingerprints import compute_docked_pose_interactions
from syndesis.stage5.final_pose_labeling import label_final_poses
from syndesis.stage5.interaction_recovery import compute_interaction_recovery
from syndesis.stage5.load_stage_inputs import load_stage5_config, load_stage5_inputs, stage5_paths
from syndesis.stage5.native_atlas import build_native_interaction_atlas as build_native_interaction_atlas_impl
from syndesis.stage5.plip_engine import run_plip_crosscheck as run_plip_crosscheck_table
from syndesis.stage5.report_stage5 import write_stage5_report
from syndesis.stage5.residue_mapping import build_interaction_residue_map


def _load(config_path: str | Path) -> tuple[dict[str, Any], dict[str, Path], dict[str, pd.DataFrame]]:
    config = load_stage5_config(config_path)
    paths = stage5_paths(config)
    inputs = load_stage5_inputs(config)
    return config, paths, inputs


def build_native_interaction_atlas(config_path: str | Path) -> dict[str, Any]:
    config, paths, inputs = _load(config_path)
    residue_map = build_interaction_residue_map(inputs, paths)
    complexes = build_interaction_complexes(inputs, paths)
    long, fps, key = build_native_interaction_atlas_impl(complexes, residue_map, inputs, config, paths)
    return {"status": "complete", "native_interactions": int(len(long)), "native_fingerprints": int(len(fps)), "key_interactions": int(len(key))}


def compute_pose_interactions(config_path: str | Path) -> dict[str, Any]:
    config, paths, inputs = _load(config_path)
    residue_map = pd.read_parquet(paths["processed"] / "interaction_residue_map.parquet") if (paths["processed"] / "interaction_residue_map.parquet").exists() else build_interaction_residue_map(inputs, paths)
    complexes = pd.read_parquet(paths["processed"] / "interaction_complexes.parquet") if (paths["processed"] / "interaction_complexes.parquet").exists() else build_interaction_complexes(inputs, paths)
    key = pd.read_parquet(paths["processed"] / "key_egfr_interactions.parquet")
    long, fps = compute_docked_pose_interactions(complexes, residue_map, inputs, key, config, paths)
    return {"status": "complete", "pose_interactions": int(len(long)), "pose_fingerprints": int(len(fps))}


def compute_recovery(config_path: str | Path) -> dict[str, Any]:
    config, paths, inputs = _load(config_path)
    docked_fps = pd.read_parquet(paths["processed"] / "docked_pose_fingerprints.parquet")
    native_fps = pd.read_parquet(paths["processed"] / "native_interaction_fingerprints.parquet")
    key = pd.read_parquet(paths["processed"] / "key_egfr_interactions.parquet")
    recovery = compute_interaction_recovery(docked_fps, native_fps, key, inputs, config, paths)
    return {"status": "complete", "recovery_rows": int(len(recovery))}


def cluster_stage5_binding_modes(config_path: str | Path) -> dict[str, Any]:
    config, paths, inputs = _load(config_path)
    native_fps = pd.read_parquet(paths["processed"] / "native_interaction_fingerprints.parquet")
    docked_fps = pd.read_parquet(paths["processed"] / "docked_pose_fingerprints.parquet")
    clusters = cluster_binding_modes(native_fps, docked_fps, inputs, config, paths)
    return {"status": "complete", "clusters": int(clusters["cluster_id"].nunique()) if not clusters.empty else 0, "rows": int(len(clusters))}


def label_stage5_final_poses(config_path: str | Path) -> dict[str, Any]:
    config, paths, inputs = _load(config_path)
    recovery = pd.read_parquet(paths["processed"] / "interaction_recovery.parquet")
    docked_fps = pd.read_parquet(paths["processed"] / "docked_pose_fingerprints.parquet")
    clusters = pd.read_parquet(paths["processed"] / "binding_mode_clusters.parquet")
    final, features = label_final_poses(recovery, docked_fps, clusters, inputs, config, paths)
    return {"status": "complete", "final_labels": int(len(final)), "stage6_features": int(len(features))}


def run_stage5_plip_crosscheck(config_path: str | Path) -> dict[str, Any]:
    config, paths, _inputs = _load(config_path)
    plip = run_plip_crosscheck_table(config, paths)
    return {"status": "complete", "plip_rows": int(len(plip))}


def report_stage5(config_path: str | Path) -> dict[str, Any]:
    config, paths, _inputs = _load(config_path)
    native_long = pd.read_parquet(paths["processed"] / "native_interactions_long.parquet")
    key = pd.read_parquet(paths["processed"] / "key_egfr_interactions.parquet")
    docked_long = pd.read_parquet(paths["processed"] / "docked_pose_interactions_long.parquet")
    recovery = pd.read_parquet(paths["processed"] / "interaction_recovery.parquet")
    clusters = pd.read_parquet(paths["processed"] / "binding_mode_clusters.parquet")
    final = pd.read_parquet(paths["processed"] / "final_pose_labels.parquet")
    plip = pd.read_parquet(paths["processed"] / "plip_crosscheck.parquet") if (paths["processed"] / "plip_crosscheck.parquet").exists() else run_plip_crosscheck_table(config, paths)
    report = write_stage5_report(native_long, key, docked_long, recovery, clusters, final, plip, paths["reports"] / "05_interaction_atlas.html")
    payload = {
        "status": "complete",
        "native_interactions": int(len(native_long)),
        "docked_interactions": int(len(docked_long)),
        "final_labels": int(len(final)),
        "report": str(report),
    }
    write_json(paths["processed"] / "stage5_summary.json", payload)
    return payload


def build_stage5_all(config_path: str | Path) -> dict[str, Any]:
    config, paths, inputs = _load(config_path)
    residue_map = build_interaction_residue_map(inputs, paths)
    complexes = build_interaction_complexes(inputs, paths)
    native_long, native_fps, key = build_native_interaction_atlas_impl(complexes, residue_map, inputs, config, paths)
    docked_long, docked_fps = compute_docked_pose_interactions(complexes, residue_map, inputs, key, config, paths)
    recovery = compute_interaction_recovery(docked_fps, native_fps, key, inputs, config, paths)
    clusters = cluster_binding_modes(native_fps, docked_fps, inputs, config, paths)
    final, _features = label_final_poses(recovery, docked_fps, clusters, inputs, config, paths)
    plip = run_plip_crosscheck_table(config, paths)
    report = write_stage5_report(native_long, key, docked_long, recovery, clusters, final, plip, paths["reports"] / "05_interaction_atlas.html")
    payload = {
        "status": "complete",
        "native_interactions": int(len(native_long)),
        "docked_interactions": int(len(docked_long)),
        "final_labels": int(len(final)),
        "report": str(report),
    }
    write_json(paths["processed"] / "stage5_summary.json", payload)
    return payload

