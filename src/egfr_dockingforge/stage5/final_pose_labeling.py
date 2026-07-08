from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pandas as pd

from egfr_dockingforge.common.io import write_table
from egfr_dockingforge.stage5.schemas import FINAL_POSE_LABEL_COLUMNS, STAGE6_INTERACTION_FEATURE_COLUMNS


def _num(value: object) -> float | None:
    parsed = pd.to_numeric(pd.Series([value]), errors="coerce").iloc[0]
    return None if pd.isna(parsed) else float(parsed)


def classify_final_pose(row: pd.Series, config: dict[str, Any]) -> tuple[str, str, str]:
    strict = float(config.get("labeling", {}).get("strict_rmsd_threshold_angstrom", 2.0))
    relaxed = float(config.get("labeling", {}).get("relaxed_rmsd_threshold_angstrom", 3.0))
    rmsd = _num(row.get("rmsd_symmetry_corrected"))
    sanity = str(row.get("sanity_status", "")).lower()
    recovery = str(row.get("interaction_recovery_label", ""))
    recall = _num(row.get("key_interaction_recall_native"))
    if recall is None:
        recall = _num(row.get("key_interaction_recall_consensus"))
    if "fail" in sanity or bool(row.get("invalid_pose_flag", False)):
        return "invalid_pose", "high", "physical sanity failed"
    if recovery == "no_reference":
        return "no_reference_pending_review", "low", "interaction reference unavailable"
    high_recovery = recovery == "high_recovery"
    moderate_recovery = recovery in {"high_recovery", "moderate_recovery"}
    poor_recovery = recovery == "poor_recovery" or (recall is not None and recall < 0.5)
    if rmsd is not None and rmsd <= strict and high_recovery:
        return "high_confidence_native_like", "high", "strict RMSD and high interaction recovery"
    if rmsd is not None and rmsd <= relaxed and poor_recovery:
        return "rmsd_good_interactions_poor", "medium", "RMSD acceptable but key interactions poorly recovered"
    if rmsd is not None and rmsd > relaxed and moderate_recovery:
        return "rmsd_poor_interactions_good", "medium", "RMSD poor but key interactions are recovered"
    if moderate_recovery and (rmsd is None or rmsd <= relaxed or row.get("ifp_tanimoto_to_consensus") is not None):
        return "plausible_binding_mode", "medium", "consensus interaction recovery supports binding mode"
    if poor_recovery:
        return "wrong_binding_mode", "medium", "key interaction recovery poor"
    return "failed_analysis", "low", "label criteria were not met"


def label_final_poses(
    recovery: pd.DataFrame,
    docked_fps: pd.DataFrame,
    clusters: pd.DataFrame,
    inputs: dict[str, pd.DataFrame],
    config: dict[str, Any],
    paths: dict[str, Path],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    score = inputs["pose_scores"].rename(columns={"original_pose_rank": "pose_rank", "original_docking_score": "original_docking_score"})
    labels = inputs["labels"][["pose_id", "invalid_pose_flag"]]
    table = recovery.merge(labels, on="pose_id", how="left")
    rows: list[dict[str, Any]] = []
    for _, row in table.iterrows():
        final_label, confidence, reason = classify_final_pose(row, config)
        rows.append(
            {
                **{col: row.get(col) for col in FINAL_POSE_LABEL_COLUMNS if col not in {"final_pose_label", "final_pose_label_confidence", "final_pose_label_reason", "warnings_json"}},
                "final_pose_label": final_label,
                "final_pose_label_confidence": confidence,
                "final_pose_label_reason": reason,
                "warnings_json": row.get("warnings_json", "[]"),
            }
        )
    final = pd.DataFrame(rows, columns=FINAL_POSE_LABEL_COLUMNS)
    write_table(paths["processed"] / "final_pose_labels.parquet", final)
    write_table(paths["processed"] / "final_pose_labels.csv", final)

    fp_cols = ["pose_id", "num_interactions", "num_key_interactions", "fingerprint_sparse_json"]
    cluster_map = clusters[clusters["entity_type"].eq("docked_pose")][["pose_id", "cluster_id"]].rename(columns={"cluster_id": "binding_mode_cluster_id"})
    features = (
        score.merge(final[["pose_id", "final_pose_label"]], on="pose_id", how="left")
        .merge(recovery[["pose_id", "ifp_tanimoto_to_native", "ifp_tanimoto_to_consensus", "key_interaction_recall_native", "key_interaction_precision_native", "key_interaction_f1_native", "key_interaction_recall_consensus", "key_interaction_precision_consensus", "key_interaction_f1_consensus"]], on="pose_id", how="left")
        .merge(docked_fps[fp_cols], on="pose_id", how="left")
        .merge(cluster_map, on="pose_id", how="left")
    )
    features = features.rename(columns={"cnnscore": "gnina_cnnscore", "cnnaffinity": "gnina_cnnaffinity"})
    for col in STAGE6_INTERACTION_FEATURE_COLUMNS:
        if col not in features.columns:
            features[col] = None
    features = features[STAGE6_INTERACTION_FEATURE_COLUMNS]
    write_table(paths["processed"] / "stage6_pose_features_interactions.parquet", features)
    write_table(paths["processed"] / "stage6_pose_features_interactions.csv", features)
    return final, features
