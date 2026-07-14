from __future__ import annotations

import pandas as pd

from syndesis.stage6.feature_builder import training_feature_columns


def test_training_feature_columns_excludes_metadata_and_json() -> None:
    features = pd.DataFrame({"pose_id": ["p1"], "cnnscore": [0.8], "fingerprint_sparse_json": ["[]"]})
    audit = pd.DataFrame({"feature_name": ["pose_id", "cnnscore", "fingerprint_sparse_json"], "allowed_for_training": [False, True, False]})
    assert training_feature_columns(features, audit) == ["cnnscore"]
