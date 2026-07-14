from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from syndesis.common.io import write_table
from syndesis.stage4.schemas import FAILURE_COLUMNS, POSE_SCORE_COLUMNS, SCORER_SUMMARY_COLUMNS, STAGE6_FEATURE_COLUMNS, TASK_METRIC_COLUMNS
from syndesis.stage4.score_normalization import sort_by_scorer


SCORERS = {
    "original_score": "original_docking_score",
    "cnnscore": "cnnscore",
    "cnnaffinity": "cnnaffinity",
    "gnina_empirical": "gnina_empirical_affinity",
}


def _bool_series(frame: pd.DataFrame, column: str) -> pd.Series:
    return frame[column].fillna(False).astype(bool) if column in frame else pd.Series(False, index=frame.index)


def build_pose_score_table(tasks: pd.DataFrame, labels: pd.DataFrame, empirical: pd.DataFrame, gnina: pd.DataFrame, config: dict[str, Any], paths: dict[str, Path]) -> pd.DataFrame:
    label_cols = [
        "pose_id", "strict_native_like_flag", "relaxed_native_like_flag", "invalid_pose_flag",
        "sampled_not_ranked_flag", "ranking_failure_flag", "sampling_failure_flag",
        "interaction_recovery_status", "final_pose_label_status",
    ]
    table = tasks.merge(labels[label_cols], on="pose_id", how="left")
    table = table.merge(empirical, on=["pose_id", "original_docking_score"], how="left")
    table = table.merge(gnina[["pose_id", "cnnscore", "cnnaffinity", "cnn_vs", "rescoring_status", "rescoring_warnings_json"]], on="pose_id", how="left")
    table["stage3_pose_label"] = table["native_like_label_stage3"]
    if "receptor_state" not in table:
        table["receptor_state"] = None
    table["interaction_recovery_status"] = "pending_stage5"
    table["final_pose_label_status"] = "pending_stage5"
    table["warnings_json"] = table["rescoring_warnings_json"].fillna("[]")
    table = table.rename(columns={"rescoring_status": "rescoring_status"})
    for col in POSE_SCORE_COLUMNS:
        if col not in table.columns:
            table[col] = None
    table = table[POSE_SCORE_COLUMNS]
    write_table(paths["processed"] / "pose_score_table.parquet", table)
    write_table(paths["processed"] / "pose_score_table.csv", table)
    features = table[STAGE6_FEATURE_COLUMNS].copy()
    write_table(paths["processed"] / "stage6_pose_features_base.parquet", features)
    write_table(paths["processed"] / "stage6_pose_features_base.csv", features)
    return table


def _rank_of_first(sorted_frame: pd.DataFrame, flag_column: str) -> int | None:
    flags = _bool_series(sorted_frame, flag_column).tolist()
    for idx, value in enumerate(flags, start=1):
        if value:
            return idx
    return None


def _auc_proxy(sorted_frame: pd.DataFrame, flag_column: str) -> float | None:
    flags = _bool_series(sorted_frame, flag_column).astype(int).to_numpy()
    if len(flags) == 0 or flags.sum() == 0 or flags.sum() == len(flags):
        return None
    ranks = np.arange(1, len(flags) + 1)
    pos_ranks = ranks[flags == 1]
    n_pos = int(flags.sum())
    n_neg = int(len(flags) - n_pos)
    auc = (n_pos * n_neg + n_pos * (n_pos + 1) / 2 - pos_ranks.sum()) / (n_pos * n_neg)
    return float(auc)


def _spearman(frame: pd.DataFrame, scorer: str, directions: dict[str, Any]) -> float | None:
    valid = frame[[scorer, "rmsd_symmetry_corrected"]].apply(pd.to_numeric, errors="coerce").dropna()
    if len(valid) < 3:
        return None
    score = valid[scorer]
    if str(directions.get(scorer, "lower")).lower() == "higher":
        score = -score
    return float(score.corr(valid["rmsd_symmetry_corrected"], method="spearman"))


def rescoring_task_metrics(score_table: pd.DataFrame, config: dict[str, Any], paths: dict[str, Path]) -> pd.DataFrame:
    directions = config["ranking"]["directions"]
    rows: list[dict[str, Any]] = []
    for task_id, group in score_table.groupby("docking_task_id", dropna=False):
        row: dict[str, Any] = {
            "docking_task_id": task_id,
            "task_type": group["task_type"].iloc[0],
            "ligand_id": group["ligand_id"].iloc[0],
            "target_receptor_id": group["target_receptor_id"].iloc[0],
            "docking_engine": group["docking_engine"].iloc[0],
            "num_poses_scored": int(group["rescoring_status"].eq("success").sum()),
            "warnings_json": json.dumps([]),
        }
        best_rmsd_by_scorer: dict[str, float | None] = {}
        for prefix, scorer in SCORERS.items():
            ordered = sort_by_scorer(group, scorer, directions)
            top = ordered.iloc[0] if not ordered.empty else None
            top_rmsd = None if top is None else pd.to_numeric(pd.Series([top["rmsd_symmetry_corrected"]]), errors="coerce").iloc[0]
            best_rmsd_by_scorer[prefix] = None if pd.isna(top_rmsd) else float(top_rmsd)
        row["top1_rmsd_original_score"] = best_rmsd_by_scorer["original_score"]
        row["top1_rmsd_cnnscore"] = best_rmsd_by_scorer["cnnscore"]
        row["top1_rmsd_cnnaffinity"] = best_rmsd_by_scorer["cnnaffinity"]
        row["top1_rmsd_gnina_empirical"] = best_rmsd_by_scorer["gnina_empirical"]
        for n in [5, 10]:
            for name, scorer in [("original_score", "original_docking_score"), ("cnnscore", "cnnscore"), ("cnnaffinity", "cnnaffinity")]:
                ordered = sort_by_scorer(group, scorer, directions).head(n)
                rmsd = pd.to_numeric(ordered["rmsd_symmetry_corrected"], errors="coerce")
                row[f"best_rmsd_top{n}_{name}"] = None if rmsd.dropna().empty else float(rmsd.min())
        strict = float(config["ranking"]["strict_rmsd_threshold_angstrom"])
        relaxed = float(config["ranking"]["relaxed_rmsd_threshold_angstrom"])
        for name, col in [("original_score", "top1_rmsd_original_score"), ("cnnscore", "top1_rmsd_cnnscore"), ("cnnaffinity", "top1_rmsd_cnnaffinity")]:
            value = row[col]
            row[f"top1_success_{name}_strict"] = bool(value is not None and value <= strict)
            row[f"top1_success_{name}_relaxed"] = bool(value is not None and value <= relaxed)
        for name, scorer in [("original_score", "original_docking_score"), ("cnnscore", "cnnscore"), ("cnnaffinity", "cnnaffinity")]:
            ordered = sort_by_scorer(group, scorer, directions)
            row[f"rank_of_first_strict_native_like_by_{name}"] = _rank_of_first(ordered, "strict_native_like_flag")
            row[f"score_native_like_auc_{name.replace('_score','')}"] = _auc_proxy(ordered, "strict_native_like_flag")
            row[f"spearman_score_vs_rmsd_{name.replace('_score','')}"] = _spearman(group, scorer, directions)
        candidates = {k: v for k, v in best_rmsd_by_scorer.items() if v is not None}
        row["dominant_best_scorer"] = min(candidates, key=candidates.get) if candidates else None
        rows.append(row)
    frame = pd.DataFrame(rows, columns=TASK_METRIC_COLUMNS)
    write_table(paths["processed"] / "rescoring_task_metrics.parquet", frame)
    write_table(paths["processed"] / "rescoring_task_metrics.csv", frame)
    return frame


def scorer_summary(score_table: pd.DataFrame, metrics: pd.DataFrame, raw_runs: pd.DataFrame, config: dict[str, Any], paths: dict[str, Path]) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    scorer_map = {
        "original_docking_score": "original_score",
        "cnnscore": "cnnscore",
        "cnnaffinity": "cnnaffinity",
    }
    for scorer, suffix in scorer_map.items():
        task_col_strict = f"top1_success_{suffix}_strict"
        task_col_relaxed = f"top1_success_{suffix}_relaxed"
        rank_col = f"rank_of_first_strict_native_like_by_{suffix}"
        auc_col = f"score_native_like_auc_{suffix.replace('_score','')}"
        rho_col = f"spearman_score_vs_rmsd_{suffix.replace('_score','')}"
        missing = pd.to_numeric(score_table[scorer], errors="coerce").isna() if scorer in score_table else pd.Series(True, index=score_table.index)
        row = {
            "scorer_name": scorer,
            "num_tasks_evaluated": int(len(metrics)),
            "redocking_top1_success_strict": _mean_for_task(metrics, task_col_strict, "redocking"),
            "redocking_top1_success_relaxed": _mean_for_task(metrics, task_col_relaxed, "redocking"),
            "same_state_crossdock_top1_success_strict": _mean_for_task(metrics, task_col_strict, "crossdocking_same_state"),
            "same_state_crossdock_top1_success_relaxed": _mean_for_task(metrics, task_col_relaxed, "crossdocking_same_state"),
            "other_state_crossdock_top1_success_strict": _mean_for_task(metrics, task_col_strict, "crossdocking_other_state"),
            "other_state_crossdock_top1_success_relaxed": _mean_for_task(metrics, task_col_relaxed, "crossdocking_other_state"),
            "mean_rank_of_first_strict_native_like": pd.to_numeric(metrics.get(rank_col), errors="coerce").mean(),
            "median_rank_of_first_strict_native_like": pd.to_numeric(metrics.get(rank_col), errors="coerce").median(),
            "score_native_like_auc_mean": pd.to_numeric(metrics.get(auc_col), errors="coerce").mean(),
            "score_native_like_auc_median": pd.to_numeric(metrics.get(auc_col), errors="coerce").median(),
            "spearman_score_vs_rmsd_median": pd.to_numeric(metrics.get(rho_col), errors="coerce").median(),
            "invalid_pose_preference_rate": 0.0,
            "missing_score_rate": float(missing.mean()) if len(missing) else None,
            "runtime_seconds_total": float(pd.to_numeric(raw_runs.get("runtime_seconds"), errors="coerce").sum()) if scorer != "original_docking_score" else 0.0,
            "runtime_seconds_per_pose_median": float(pd.to_numeric(raw_runs.get("runtime_seconds"), errors="coerce").median()) if scorer != "original_docking_score" and not raw_runs.empty else 0.0,
        }
        rows.append(row)
    frame = pd.DataFrame(rows, columns=SCORER_SUMMARY_COLUMNS)
    write_table(paths["processed"] / "scorer_comparison_summary.parquet", frame)
    write_table(paths["processed"] / "scorer_comparison_summary.csv", frame)
    return frame


def _mean_for_task(metrics: pd.DataFrame, column: str, task_pattern: str) -> float | None:
    if column not in metrics:
        return None
    subset = metrics[metrics["task_type"].astype(str).str.contains(task_pattern, regex=False)]
    if subset.empty:
        return None
    return float(subset[column].fillna(False).astype(bool).mean())


def failure_table(score_table: pd.DataFrame, paths: dict[str, Path]) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for task_id, group in score_table.groupby("docking_task_id", dropna=False):
        if group["rescoring_status"].ne("success").any():
            failed = group[group["rescoring_status"].ne("success")]
            for _, row in failed.iterrows():
                rows.append(_failure(row, "gnina_failed_to_score", "GNINA failed or produced no successful score.", "Check logs and input files."))
        for scorer, failure_type in [("cnnscore", "cnnscore_prefers_wrong_pose"), ("cnnaffinity", "cnnaffinity_prefers_wrong_pose")]:
            ordered = group.assign(_score=pd.to_numeric(group[scorer], errors="coerce")).sort_values("_score", ascending=False, na_position="last")
            if not ordered.empty and not bool(ordered.iloc[0].get("relaxed_native_like_flag", False)):
                rows.append(_failure(ordered.iloc[0], failure_type, f"{scorer} top-ranked a non-native-like Stage 3 pose.", "Inspect Stage 5 interactions and receptor state."))
    frame = pd.DataFrame(rows, columns=FAILURE_COLUMNS)
    write_table(paths["processed"] / "rescoring_failures.parquet", frame)
    write_table(paths["processed"] / "rescoring_failures.csv", frame)
    return frame


def _failure(row: pd.Series, failure_type: str, description: str, followup: str) -> dict[str, Any]:
    return {
        "pose_id": row.get("pose_id"),
        "docking_task_id": row.get("docking_task_id"),
        "ligand_id": row.get("ligand_id"),
        "target_receptor_id": row.get("target_receptor_id"),
        "task_type": row.get("task_type"),
        "failure_type": failure_type,
        "failure_description": description,
        "original_docking_score": row.get("original_docking_score"),
        "cnnscore": row.get("cnnscore"),
        "cnnaffinity": row.get("cnnaffinity"),
        "rmsd_symmetry_corrected": row.get("rmsd_symmetry_corrected"),
        "sanity_status": row.get("sanity_status"),
        "stage3_pose_label": row.get("stage3_pose_label"),
        "suggested_followup": followup,
    }
