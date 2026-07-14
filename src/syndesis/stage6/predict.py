from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import joblib
import numpy as np
import pandas as pd
from pandas.api.types import is_string_dtype


def _prepare_x(feature_table: pd.DataFrame, feature_columns: list[str], categorical_columns: list[str]) -> pd.DataFrame:
    x = feature_table.reindex(columns=feature_columns).copy()
    for col in categorical_columns:
        if col in x:
            x[col] = x[col].astype("category").cat.codes
    for col in x.columns:
        if x[col].dtype == object or is_string_dtype(x[col]):
            x[col] = x[col].astype("category").cat.codes
    return x.replace([np.inf, -np.inf], np.nan).fillna(0)


def score_poses(feature_table: pd.DataFrame, model_artifact: str | Path | dict[str, Any]) -> pd.DataFrame:
    artifact = joblib.load(model_artifact) if not isinstance(model_artifact, dict) else model_artifact
    feature_columns = artifact["feature_columns"]
    categorical_columns = artifact.get("categorical_columns", [])
    x = _prepare_x(feature_table, feature_columns, categorical_columns)
    ranker = artifact["ranker"]
    classifier = artifact.get("calibrator") or artifact.get("confidence_classifier")
    rank_score = ranker.predict(x)
    if classifier is not None:
        prob = classifier.predict_proba(x)[:, 1]
    else:
        prob = np.full(len(x), np.nan)
    out = feature_table[["pose_id", "docking_task_id"]].rename(columns={"docking_task_id": "group_id"}).copy()
    out["pose_rank_model_score"] = rank_score
    out["pose_confidence_probability"] = prob
    out["model_rank_within_group"] = out.groupby("group_id")["pose_rank_model_score"].rank(ascending=False, method="first").astype(int)
    out["confidence_label"] = pd.cut(
        out["pose_confidence_probability"],
        bins=[-0.01, 0.3, 0.7, 1.01],
        labels=["low", "medium", "high"],
    ).astype(str)
    out["explanation_summary_json"] = [
        json.dumps({"top_driver": "see_stage6_feature_importance", "rank_score": float(score)})
        for score in out["pose_rank_model_score"]
    ]
    return out
