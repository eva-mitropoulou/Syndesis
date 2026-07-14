from __future__ import annotations

import pandas as pd

from syndesis.stage6.model_selection import prepare_model_matrix
from syndesis.stage6.xgb_ranker import make_xgb_ranker


def test_prepare_model_matrix_encodes_categoricals() -> None:
    x, categoricals = prepare_model_matrix(pd.DataFrame({"a": ["x", "y"], "b": [1.0, None]}), ["a", "b"])
    assert categoricals == ["a"]
    assert x.isna().sum().sum() == 0
    assert x["a"].dtype.kind in "iu"


def test_xgb_ranker_uses_graded_relevance_objective() -> None:
    model = make_xgb_ranker("rank:ndcg", seed=1, n_jobs=1, n_estimators=2)
    assert model.get_params()["objective"] == "rank:ndcg"
