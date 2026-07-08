from __future__ import annotations

import json
import platform
import time
from pathlib import Path
from typing import Any

import importlib.util
import joblib
import numpy as np
import pandas as pd
from pandas.api.types import is_string_dtype
from sklearn.metrics import (
    average_precision_score,
    brier_score_loss,
    ndcg_score,
    roc_auc_score,
)

from egfr_dockingforge.common.io import write_json, write_table
from egfr_dockingforge.stage6.baselines import baseline_scores, classifier_baselines, rule_based_probability
from egfr_dockingforge.stage6.calibration import calibrate_classifier, expected_calibration_error
from egfr_dockingforge.stage6.catboost_ranker import make_catboost_ranker
from egfr_dockingforge.stage6.confidence_classifier import make_lgbm_classifier, make_xgb_classifier
from egfr_dockingforge.stage6.feature_builder import training_feature_columns
from egfr_dockingforge.stage6.leakage_audit import assert_no_leakage
from egfr_dockingforge.stage6.lgbm_ranker import make_lgbm_ranker
from egfr_dockingforge.stage6.xgb_ranker import make_xgb_ranker


def prepare_model_matrix(features: pd.DataFrame, feature_columns: list[str]) -> tuple[pd.DataFrame, list[str]]:
    x = features.reindex(columns=feature_columns).copy()
    categorical_columns: list[str] = []
    for col in x.columns:
        if x[col].dtype == object or is_string_dtype(x[col]):
            categorical_columns.append(col)
            x[col] = x[col].astype("category").cat.codes
    x = x.replace([np.inf, -np.inf], np.nan).fillna(0)
    return x, categorical_columns


def _group_sizes(group_ids: pd.Series) -> list[int]:
    return group_ids.value_counts(sort=False).to_list()


def _ordered_subset(x: pd.DataFrame, y: pd.Series, meta: pd.DataFrame, mask: pd.Series) -> tuple[pd.DataFrame, pd.Series, pd.DataFrame]:
    idx = meta.loc[mask].sort_values(["group_id", "pose_id"]).index
    return x.loc[idx], y.loc[idx], meta.loc[idx]


def _ranking_metrics(scores: pd.Series | np.ndarray, labels: pd.Series, meta: pd.DataFrame, model_id: str, split_name: str, tvt: str, family: str) -> list[dict[str, Any]]:
    tmp = meta[["pose_id", "group_id", "receptor_state"]].copy()
    tmp["score"] = np.asarray(scores, dtype=float)
    tmp["label"] = labels.to_numpy(dtype=int)
    rows = []
    ndcg1: list[float] = []
    ndcg3: list[float] = []
    reciprocal: list[float] = []
    top1: list[float] = []
    top3: list[float] = []
    first_ranks: list[float] = []
    for _, group in tmp.groupby("group_id"):
        if group["label"].nunique() <= 1:
            continue
        y_true = group["label"].to_numpy(dtype=float).reshape(1, -1)
        y_score = group["score"].to_numpy(dtype=float).reshape(1, -1)
        ndcg1.append(float(ndcg_score(y_true, y_score, k=1)))
        ndcg3.append(float(ndcg_score(y_true, y_score, k=min(3, len(group)))))
        ranked = group.sort_values("score", ascending=False).reset_index(drop=True)
        relevant = ranked.index[ranked["label"].gt(0)].to_list()
        if relevant:
            first = relevant[0] + 1
            reciprocal.append(1.0 / first)
            first_ranks.append(float(first))
            top1.append(float(first == 1))
            top3.append(float(first <= 3))
    def add(name: str, values: list[float]) -> None:
        rows.append(
            {
                "model_id": model_id,
                "model_family": family,
                "split_name": split_name,
                "split_fold": 0,
                "metric_name": f"{tvt}_{name}",
                "metric_value": float(np.mean(values)) if values else np.nan,
                "metric_context": "ranking",
                "num_groups": int(tmp["group_id"].nunique()),
                "num_poses": int(len(tmp)),
                "notes": "",
            }
        )
    add("NDCG@1", ndcg1)
    add("NDCG@3", ndcg3)
    add("MRR", reciprocal)
    add("top1_high_confidence_pose_rate", top1)
    add("top3_contains_high_confidence_pose_rate", top3)
    add("mean_rank_of_first_high_confidence_pose", first_ranks)
    return rows


def _classifier_metrics(prob: np.ndarray, labels: pd.Series, meta: pd.DataFrame, model_id: str, split_name: str, tvt: str, family: str) -> list[dict[str, Any]]:
    y = labels.to_numpy(dtype=int)
    rows = []
    values = {
        "Brier": brier_score_loss(y, prob),
        "expected_calibration_error": expected_calibration_error(y, prob),
    }
    if len(set(y)) == 2:
        values["ROC-AUC"] = roc_auc_score(y, prob)
        values["PR-AUC"] = average_precision_score(y, prob)
    for threshold in (0.3, 0.5, 0.7):
        pred = prob >= threshold
        values[f"precision_at_{threshold}"] = float(y[pred].mean()) if pred.any() else np.nan
        values[f"recall_at_{threshold}"] = float(pred[y == 1].mean()) if (y == 1).any() else np.nan
    for name, value in values.items():
        rows.append(
            {
                "model_id": model_id,
                "model_family": family,
                "split_name": split_name,
                "split_fold": 0,
                "metric_name": f"{tvt}_{name}",
                "metric_value": float(value),
                "metric_context": "classification",
                "num_groups": int(meta["group_id"].nunique()),
                "num_poses": int(len(meta)),
                "notes": "",
            }
        )
    return rows


def _primary_metric(metrics: pd.DataFrame, model_id: str) -> float:
    row = metrics[(metrics["model_id"].eq(model_id)) & (metrics["metric_name"].eq("valid_NDCG@1"))]
    if row.empty:
        return float("-inf")
    return float(row["metric_value"].iloc[0])


def train_and_evaluate_models(features: pd.DataFrame, labels: pd.DataFrame, groups: pd.DataFrame, splits: pd.DataFrame, audit: pd.DataFrame, config: dict, paths: dict) -> dict[str, Any]:
    start = time.time()
    feature_columns = training_feature_columns(features, audit)
    assert_no_leakage(feature_columns)
    x, categorical_columns = prepare_model_matrix(features, feature_columns)
    y_rank = labels.set_index("pose_id").loc[features["pose_id"], "rank_relevance_label"].astype(int).reset_index(drop=True)
    y_bin = labels.set_index("pose_id").loc[features["pose_id"], "binary_confidence_label"].astype(int).reset_index(drop=True)
    meta = features[["pose_id", "docking_task_id", "ligand_id", "target_receptor_id", "receptor_state"]].copy()
    meta["group_id"] = meta["docking_task_id"]
    usable_groups = set(groups.loc[groups["group_usable_for_ranking"], "group_id"])
    primary_split = config["splits"]["primary_split_name"]
    split_rows = splits[splits["split_name"].eq(primary_split)].set_index("pose_id").loc[features["pose_id"]].reset_index()
    meta["train_valid_test"] = split_rows["train_valid_test"].to_numpy()

    train_mask = meta["train_valid_test"].eq("train") & meta["group_id"].isin(usable_groups)
    valid_mask = meta["train_valid_test"].eq("valid")
    test_mask = meta["train_valid_test"].eq("test")
    if y_bin[meta["train_valid_test"].eq("train")].nunique() < 2:
        raise RuntimeError("Stage 6 confidence training split has only one binary class.")
    if not train_mask.any():
        raise RuntimeError("Stage 6 ranker training has no usable groups with relevance variation.")

    metrics_rows: list[dict[str, Any]] = []
    ranker_artifacts: dict[str, Any] = {}
    model_cfg = config["models"]
    seed = int(model_cfg["seed"])
    n_jobs = int(model_cfg["n_jobs"])
    n_estimators = int(model_cfg["ranker_estimators"])

    for model_id, score in baseline_scores(features).items():
        for tvt, mask in [("valid", valid_mask), ("test", test_mask)]:
            metrics_rows.extend(_ranking_metrics(score.loc[mask], y_rank.loc[mask], meta.loc[mask], model_id, primary_split, tvt, "baseline_ranker"))

    ranker_specs = [
        ("lgbm_rank_xendcg", make_lgbm_ranker("rank_xendcg", seed, n_jobs, n_estimators), "lightgbm_ranker"),
        ("lgbm_lambdarank", make_lgbm_ranker("lambdarank", seed, n_jobs, n_estimators), "lightgbm_ranker"),
        ("xgb_rank_ndcg", make_xgb_ranker("rank:ndcg", seed, n_jobs, n_estimators), "xgboost_ranker"),
    ]
    if importlib.util.find_spec("catboost") is not None:
        ranker_specs.append(("catboost_yetirank", make_catboost_ranker("YetiRank", seed, n_estimators), "catboost_ranker"))
    x_train, y_train, meta_train = _ordered_subset(x, y_rank, meta, train_mask)
    x_valid, y_valid_rank, meta_valid = _ordered_subset(x, y_rank, meta, valid_mask)
    for model_id, model, family in ranker_specs:
        fit_start = time.time()
        if family == "lightgbm_ranker":
            model.fit(x_train, y_train, group=_group_sizes(meta_train["group_id"]), eval_set=[(x_valid, y_valid_rank)], eval_group=[_group_sizes(meta_valid["group_id"])])
        elif family == "xgboost_ranker":
            train_qid = pd.Series(pd.factorize(meta_train["group_id"])[0], index=meta_train.index)
            valid_levels = {value: idx for idx, value in enumerate(sorted(meta_valid["group_id"].unique()))}
            valid_qid = meta_valid["group_id"].map(valid_levels).astype(int)
            model.fit(x_train, y_train, qid=train_qid, eval_set=[(x_valid, y_valid_rank)], eval_qid=[valid_qid], verbose=False)
        elif family == "catboost_ranker":
            from catboost import Pool

            model.fit(Pool(x_train, y_train, group_id=meta_train["group_id"]), eval_set=Pool(x_valid, y_valid_rank, group_id=meta_valid["group_id"]))
        ranker_artifacts[model_id] = model
        for tvt, mask in [("valid", valid_mask), ("test", test_mask)]:
            metrics_rows.extend(_ranking_metrics(model.predict(x.loc[mask]), y_rank.loc[mask], meta.loc[mask], model_id, primary_split, tvt, family))
        metrics_rows.append(
            {
                "model_id": model_id,
                "model_family": family,
                "split_name": primary_split,
                "split_fold": 0,
                "metric_name": "runtime_seconds",
                "metric_value": float(time.time() - fit_start),
                "metric_context": "operational",
                "num_groups": int(meta_train["group_id"].nunique()),
                "num_poses": int(len(meta_train)),
                "notes": "",
            }
        )

    classifier_artifacts: dict[str, Any] = {}
    train_cls = meta["train_valid_test"].eq("train")
    classifier_specs = classifier_baselines(seed, n_jobs)
    classifier_specs["lgbm_binary_classifier"] = make_lgbm_classifier(seed, n_jobs, int(model_cfg["classifier_estimators"]))
    classifier_specs["xgb_binary_classifier"] = make_xgb_classifier(seed, n_jobs, int(model_cfg["classifier_estimators"]))
    x_train_cls, y_train_cls, _ = _ordered_subset(x, y_bin, meta, train_cls)
    x_valid_cls, y_valid_cls, meta_valid_cls = _ordered_subset(x, y_bin, meta, valid_mask)
    for model_id, model in classifier_specs.items():
        model.fit(x_train_cls, y_train_cls)
        classifier_artifacts[model_id] = model
        for tvt, mask in [("valid", valid_mask), ("test", test_mask)]:
            prob = model.predict_proba(x.loc[mask])[:, 1]
            metrics_rows.extend(_classifier_metrics(prob, y_bin.loc[mask], meta.loc[mask], model_id, primary_split, tvt, "confidence_classifier"))
    rule_prob = rule_based_probability(features)
    for tvt, mask in [("valid", valid_mask), ("test", test_mask)]:
        metrics_rows.extend(_classifier_metrics(rule_prob[mask.to_numpy()], y_bin.loc[mask], meta.loc[mask], "simple_rule_based_threshold_model", primary_split, tvt, "baseline_classifier"))

    metrics = pd.DataFrame(metrics_rows)
    write_table(paths["processed"] / "model_metrics.parquet", metrics)
    write_table(paths["processed"] / "model_metrics.csv", metrics)

    ranker_scores = {model_id: _primary_metric(metrics, model_id) for model_id in ranker_artifacts}
    best_ranker_id = max(ranker_scores, key=ranker_scores.get)
    confidence_scores = {
        model_id: float(metrics[(metrics["model_id"].eq(model_id)) & (metrics["metric_name"].eq("valid_PR-AUC"))]["metric_value"].max())
        for model_id in classifier_artifacts
    }
    best_classifier_id = max(confidence_scores, key=confidence_scores.get)
    calibrated, calibration_metrics = calibrate_classifier(classifier_artifacts[best_classifier_id], x_valid_cls, y_valid_cls, best_classifier_id, paths)

    artifact = {
        "ranker": ranker_artifacts[best_ranker_id],
        "confidence_classifier": classifier_artifacts[best_classifier_id],
        "calibrator": calibrated,
        "feature_columns": feature_columns,
        "categorical_columns": categorical_columns,
        "selected_ranker_model_id": best_ranker_id,
        "selected_confidence_model_id": best_classifier_id,
        "library_versions": _library_versions(),
    }
    artifact_path = paths["models"] / "pose_ranker.pkl"
    joblib.dump(artifact, artifact_path)
    joblib.dump(classifier_artifacts[best_classifier_id], paths["models"] / "pose_confidence_classifier.pkl")
    joblib.dump(calibrated, paths["models"] / "pose_confidence_calibrator.pkl")
    write_json(paths["models"] / "feature_schema.json", {"feature_columns": feature_columns, "categorical_columns": categorical_columns})

    docking = _primary_metric(metrics, "original_docking_score_ranker")
    gnina = _primary_metric(metrics, "gnina_cnnscore_ranker")
    best_value = ranker_scores[best_ranker_id]
    passed = best_value > docking and best_value > gnina
    selection_reason = "selected_model_beats_docking_and_gnina_baselines" if passed else "selection_failed_best_model_did_not_beat_required_baselines"
    selection = pd.DataFrame(
        [
            {
                "selected_model_id": "stage6_pose_ranker_confidence_v1",
                "selected_ranker_model_id": best_ranker_id,
                "selected_confidence_model_id": best_classifier_id,
                "selection_reason": selection_reason,
                "primary_validation_metric": best_value,
                "baseline_comparison_summary": json.dumps({"docking": docking, "gnina_cnnscore": gnina, "best": best_value}),
                "known_limitations": "Small EGFR cocrystal-derived benchmark; use Stage 7/9 review before prospective decisions.",
                "model_artifact_path": str(artifact_path),
                "calibration_artifact_path": str(paths["models"] / "pose_confidence_calibrator.pkl"),
                "feature_list_path": str(paths["models"] / "feature_schema.json"),
            }
        ]
    )
    write_table(paths["processed"] / "model_selection_summary.parquet", selection)
    write_table(paths["processed"] / "model_selection_summary.csv", selection)
    _write_model_card(paths["models"] / "model_card.md", selection.iloc[0].to_dict(), len(features), len(feature_columns), time.time() - start)
    write_json(paths["processed"] / "stage6_training_summary.json", {"status": "complete", "selected_ranker": best_ranker_id, "selected_classifier": best_classifier_id, "selection_passed": passed})
    return {
        "metrics": metrics,
        "selection": selection,
        "calibration_metrics": calibration_metrics,
        "artifact_path": artifact_path,
        "feature_columns": feature_columns,
        "categorical_columns": categorical_columns,
    }


def _library_versions() -> dict[str, str]:
    import lightgbm
    import sklearn
    import xgboost

    versions = {
        "python": platform.python_version(),
        "lightgbm": lightgbm.__version__,
        "xgboost": xgboost.__version__,
        "scikit_learn": sklearn.__version__,
    }
    if importlib.util.find_spec("catboost") is not None:
        import catboost

        versions["catboost"] = catboost.__version__
    return versions


def _write_model_card(path: Path, selection: dict[str, Any], rows: int, features: int, runtime: float) -> None:
    text = f"""# Stage 6 Pose Reranking Confidence Model

Rows: {rows}
Training features: {features}
Selected ranker: {selection['selected_ranker_model_id']}
Selected confidence model: {selection['selected_confidence_model_id']}
Selection reason: {selection['selection_reason']}
Primary validation metric: {selection['primary_validation_metric']}
Runtime seconds: {runtime:.2f}

This artifact ranks poses within docking tasks and produces calibrated pose-confidence probabilities. It is not an activity model, MD workflow, analog generator, or vendor screen.
"""
    path.write_text(text, encoding="utf-8")
