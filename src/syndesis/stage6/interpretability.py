from __future__ import annotations

import json

import joblib
import numpy as np
import pandas as pd
import shap
from sklearn.inspection import permutation_importance
from sklearn.metrics import ndcg_score

from syndesis.common.io import write_table
from syndesis.stage6.model_selection import prepare_model_matrix


def write_interpretability(features: pd.DataFrame, labels: pd.DataFrame, audit: pd.DataFrame, paths: dict) -> tuple[pd.DataFrame, pd.DataFrame]:
    artifact = joblib.load(paths["models"] / "pose_ranker.pkl")
    feature_columns = artifact["feature_columns"]
    x, _ = prepare_model_matrix(features, feature_columns)
    y = labels.set_index("pose_id").loc[features["pose_id"], "rank_relevance_label"].astype(int)
    ranker = artifact["ranker"]
    importances = getattr(ranker, "feature_importances_", np.zeros(len(feature_columns)))
    fi = pd.DataFrame({"feature_name": feature_columns, "importance": importances})
    fi = fi.merge(audit[["feature_name", "feature_group"]], on="feature_name", how="left").sort_values("importance", ascending=False)

    def rank_scorer(estimator, matrix, target) -> float:
        if len(set(target)) <= 1:
            return 0.0
        return float(ndcg_score([target], [estimator.predict(matrix)]))

    perm = permutation_importance(ranker, x, y, n_repeats=5, random_state=17, scoring=rank_scorer)
    perm_df = pd.DataFrame({"feature_name": feature_columns, "permutation_importance": perm.importances_mean})
    fi = fi.merge(perm_df, on="feature_name", how="left")
    sample = x.head(min(100, len(x)))
    explainer = shap.TreeExplainer(ranker)
    shap_values = explainer.shap_values(sample)
    shap_abs = np.abs(np.asarray(shap_values)).mean(axis=0)
    shap_df = pd.DataFrame({"feature_name": sample.columns, "mean_abs_shap": shap_abs})
    fi = fi.merge(shap_df, on="feature_name", how="left")
    write_table(paths["processed"] / "feature_importance.parquet", fi)
    write_table(paths["processed"] / "feature_importance.csv", fi)

    groups = audit[audit["allowed_for_training"]].groupby("feature_group")["feature_name"].apply(list).to_dict()
    rows = []
    full_score = float(np.nanmean(ranker.predict(x)))
    for group, cols in groups.items():
        usable = [c for c in cols if c in x.columns]
        ablated = x.copy()
        ablated[usable] = 0
        rows.append(
            {
                "feature_group": group,
                "num_features_ablated": len(usable),
                "mean_score_full": full_score,
                "mean_score_ablated": float(np.nanmean(ranker.predict(ablated))),
                "ablation_notes": json.dumps({"columns": usable[:25]}),
            }
        )
    ablation = pd.DataFrame(rows)
    write_table(paths["processed"] / "feature_group_ablation.parquet", ablation)
    write_table(paths["processed"] / "feature_group_ablation.csv", ablation)
    return fi, ablation
