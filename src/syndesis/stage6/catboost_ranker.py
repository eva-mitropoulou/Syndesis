from __future__ import annotations

def make_catboost_ranker(loss_function: str, seed: int, n_estimators: int):
    from catboost import CatBoostRanker

    return CatBoostRanker(
        loss_function=loss_function,
        iterations=n_estimators,
        learning_rate=0.05,
        depth=4,
        random_seed=seed,
        verbose=False,
        allow_writing_files=False,
    )
