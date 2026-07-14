from __future__ import annotations

import numpy as np
import pandas as pd
from sklearn.ensemble import ExtraTreesClassifier, RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler


def baseline_scores(features: pd.DataFrame) -> dict[str, pd.Series]:
    scores: dict[str, pd.Series] = {}
    if "original_docking_score" in features:
        scores["original_docking_score_ranker"] = -pd.to_numeric(features["original_docking_score"], errors="coerce").fillna(0)
    if "cnnscore" in features:
        scores["gnina_cnnscore_ranker"] = pd.to_numeric(features["cnnscore"], errors="coerce").fillna(0)
    if "cnnaffinity" in features:
        scores["gnina_cnnaffinity_ranker"] = pd.to_numeric(features["cnnaffinity"], errors="coerce").fillna(0)
    if "key_interaction_f1_consensus" in features:
        scores["prolif_consensus_ranker"] = pd.to_numeric(features["key_interaction_f1_consensus"], errors="coerce").fillna(0)
    parts = [value.rank(pct=True) for value in scores.values()]
    if parts:
        scores["simple_weighted_consensus"] = sum(parts) / len(parts)
    return scores


def classifier_baselines(seed: int, n_jobs: int) -> dict[str, object]:
    return {
        "logistic_regression_classifier": Pipeline(
            [("scale", StandardScaler()), ("model", LogisticRegression(max_iter=2000, class_weight="balanced", random_state=seed))]
        ),
        "random_forest_classifier": RandomForestClassifier(
            n_estimators=300, random_state=seed, n_jobs=n_jobs, class_weight="balanced_subsample", min_samples_leaf=2
        ),
        "extratrees_classifier": ExtraTreesClassifier(
            n_estimators=300, random_state=seed, n_jobs=n_jobs, class_weight="balanced", min_samples_leaf=2
        ),
    }


def rule_based_probability(features: pd.DataFrame) -> np.ndarray:
    cnn = pd.to_numeric(features.get("cnnscore", 0), errors="coerce").fillna(0)
    prolif = pd.to_numeric(features.get("key_interaction_f1_consensus", 0), errors="coerce").fillna(0)
    sanity = pd.to_numeric(features.get("sanity_status_encoded", 0), errors="coerce").fillna(0)
    raw = 0.45 * cnn.rank(pct=True) + 0.40 * prolif.rank(pct=True) + 0.15 * sanity
    return np.clip(raw.to_numpy(dtype=float), 0, 1)
