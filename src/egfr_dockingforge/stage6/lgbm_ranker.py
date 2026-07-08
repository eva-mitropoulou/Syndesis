from __future__ import annotations

import lightgbm as lgb


def make_lgbm_ranker(objective: str, seed: int, n_jobs: int, n_estimators: int) -> lgb.LGBMRanker:
    return lgb.LGBMRanker(
        objective=objective,
        metric="ndcg",
        n_estimators=n_estimators,
        learning_rate=0.05,
        num_leaves=15,
        min_child_samples=2,
        subsample=0.9,
        colsample_bytree=0.9,
        random_state=seed,
        n_jobs=n_jobs,
        verbose=-1,
    )
