from __future__ import annotations

from pathlib import Path
from typing import Any

import pandas as pd

from syndesis.common.io import ensure_dir, write_json, write_table
from syndesis.stage2.ensemble_export import export_ensemble
from syndesis.stage2.load_stage1 import (
    filter_stage1_candidates,
    load_stage1_benchmark,
    load_stage2_config,
    stage2_paths,
)
from syndesis.stage2.medoid_selection import select_ensemble
from syndesis.stage2.pocket_alignment import alignment_metrics, copy_aligned_receptors, update_features_with_alignment
from syndesis.stage2.pocket_clustering import cluster_receptors, distance_matrix
from syndesis.stage2.pocket_features import feature_record, features_frame
from syndesis.stage2.pocket_mapping import map_pocket_residues, mapping_frame
from syndesis.stage2.report_stage2 import write_stage2_report


def with_receptor_ids(candidates: pd.DataFrame) -> pd.DataFrame:
    frame = candidates.copy()
    frame["receptor_id"] = frame["complex_id"].astype(str).str.lower()
    frame["receptor_file_path"] = frame["receptor_clean_path"]
    return frame


def normalize_bool_columns(frame: pd.DataFrame) -> pd.DataFrame:
    frame = frame.copy()
    for column in frame.columns:
        if column.endswith("_flag") or column in {"selected_flag", "cluster_medoid_flag"}:
            try:
                frame[column] = frame[column].astype("boolean")
            except Exception:
                pass
    return frame


def build_receptor_features(config_path: str | Path) -> dict[str, Any]:
    config = load_stage2_config(config_path)
    paths = stage2_paths(config)
    for path in paths.values():
        ensure_dir(path)

    benchmark = load_stage1_benchmark(config)
    candidates = with_receptor_ids(
        filter_stage1_candidates(benchmark, config, paths["processed"] / "receptor_exclusion_table.parquet")
    )
    mapping_rows = []
    for _, row in candidates.iterrows():
        mapping_rows.extend(map_pocket_residues(row, config))
    mapping = mapping_frame(mapping_rows)
    records = [feature_record(row, config) for _, row in candidates.iterrows()]
    features = features_frame(records)

    metrics = alignment_metrics(features, config)
    features = update_features_with_alignment(features, metrics)
    aligned_lookup = copy_aligned_receptors(features, paths["aligned_receptors"], config)

    write_table(paths["processed"] / "pocket_residue_mapping.csv", mapping)
    write_table(paths["processed"] / "pocket_residue_mapping.parquet", normalize_bool_columns(mapping))
    write_table(paths["processed"] / "receptor_alignment_metrics.csv", metrics)
    write_table(paths["processed"] / "receptor_alignment_metrics.parquet", metrics)
    write_table(paths["processed"] / "aligned_receptor_lookup.csv", aligned_lookup)
    write_table(paths["processed"] / "receptor_features.csv", features)
    write_table(paths["processed"] / "receptor_features.parquet", normalize_bool_columns(features))
    return {"status": "complete", "receptor_count": int(len(features))}


def cluster_receptor_features(config_path: str | Path) -> dict[str, Any]:
    config = load_stage2_config(config_path)
    paths = stage2_paths(config)
    features_path = paths["processed"] / "receptor_features.parquet"
    if not features_path.exists():
        build_receptor_features(config_path)
    features = pd.read_parquet(features_path)
    distances = distance_matrix(features, config)
    clusters = cluster_receptors(features, distances, config)
    write_table(paths["processed"] / "receptor_distance_matrix.csv", distances)
    write_table(paths["processed"] / "receptor_distance_matrix.parquet", distances)
    write_table(paths["processed"] / "receptor_clusters.csv", clusters)
    write_table(paths["processed"] / "receptor_clusters.parquet", normalize_bool_columns(clusters))
    return {"status": "complete", "cluster_count": int(clusters["cluster_id"].nunique()) if not clusters.empty else 0}


def select_receptor_ensemble(config_path: str | Path) -> dict[str, Any]:
    config = load_stage2_config(config_path)
    paths = stage2_paths(config)
    if not (paths["processed"] / "receptor_clusters.parquet").exists():
        cluster_receptor_features(config_path)
    features = pd.read_parquet(paths["processed"] / "receptor_features.parquet")
    clusters = pd.read_parquet(paths["processed"] / "receptor_clusters.parquet")
    selected, holdout = select_ensemble(features, clusters, config)
    aligned_lookup = pd.read_csv(paths["processed"] / "aligned_receptor_lookup.csv")
    ensemble = export_ensemble(selected, aligned_lookup, paths["ensemble_receptors"], config)
    write_table(paths["processed"] / "receptor_ensemble_v1.csv", ensemble)
    write_table(paths["processed"] / "receptor_ensemble_v1.parquet", normalize_bool_columns(ensemble))
    write_table(paths["processed"] / "receptor_holdout_v1.csv", holdout)
    write_table(paths["processed"] / "receptor_holdout_v1.parquet", normalize_bool_columns(holdout))
    return {"status": "complete", "selected_receptors": int(len(ensemble)), "holdout_receptors": int(len(holdout))}


def export_receptor_ensemble(config_path: str | Path) -> dict[str, Any]:
    return select_receptor_ensemble(config_path)


def report_stage2(config_path: str | Path) -> dict[str, Any]:
    config = load_stage2_config(config_path)
    paths = stage2_paths(config)
    if not (paths["processed"] / "receptor_ensemble_v1.parquet").exists():
        select_receptor_ensemble(config_path)
    features = pd.read_parquet(paths["processed"] / "receptor_features.parquet")
    excluded = pd.read_parquet(paths["processed"] / "receptor_exclusion_table.parquet")
    clusters = pd.read_parquet(paths["processed"] / "receptor_clusters.parquet")
    ensemble = pd.read_parquet(paths["processed"] / "receptor_ensemble_v1.parquet")
    holdout = pd.read_parquet(paths["processed"] / "receptor_holdout_v1.parquet")
    report = write_stage2_report(features, excluded, clusters, ensemble, holdout, paths["reports"] / "02_receptor_ensemble.html")
    summary = {
        "status": "complete",
        "stage1_included_loaded": int(len(features) + len(excluded)),
        "stage2_passing_receptors": int(len(features)),
        "selected_receptors": int(len(ensemble)),
        "holdout_receptors": int(len(holdout)),
        "report": str(report),
    }
    write_json(paths["processed"] / "stage2_summary.json", summary)
    return summary


def build_stage2_all(config_path: str | Path) -> dict[str, Any]:
    build_receptor_features(config_path)
    cluster_receptor_features(config_path)
    select_receptor_ensemble(config_path)
    return report_stage2(config_path)
