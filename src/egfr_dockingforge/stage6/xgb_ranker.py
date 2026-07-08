from __future__ import annotations

import xgboost as xgb


def make_xgb_ranker(objective: str, seed: int, n_jobs: int, n_estimators: int) -> xgb.XGBRanker:
    return xgb.XGBRanker(
        objective=objective,
        eval_metric="ndcg@3",
        n_estimators=n_estimators,
        learning_rate=0.05,
        max_depth=4,
        min_child_weight=1,
        subsample=0.9,
        colsample_bytree=0.9,
        random_state=seed,
        n_jobs=n_jobs,
        tree_method="hist",
    )
