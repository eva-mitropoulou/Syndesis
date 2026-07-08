from __future__ import annotations

import lightgbm as lgb
import xgboost as xgb


def make_lgbm_classifier(seed: int, n_jobs: int, n_estimators: int) -> lgb.LGBMClassifier:
    return lgb.LGBMClassifier(
        objective="binary",
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


def make_xgb_classifier(seed: int, n_jobs: int, n_estimators: int) -> xgb.XGBClassifier:
    return xgb.XGBClassifier(
        objective="binary:logistic",
        eval_metric="logloss",
        n_estimators=n_estimators,
        learning_rate=0.05,
        max_depth=4,
        subsample=0.9,
        colsample_bytree=0.9,
        random_state=seed,
        n_jobs=n_jobs,
        tree_method="hist",
    )
