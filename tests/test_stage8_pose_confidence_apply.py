from __future__ import annotations

import json
from pathlib import Path


def test_stage6_feature_schema_exists() -> None:
    schema = json.loads(Path("models/stage6/feature_schema.json").read_text(encoding="utf-8"))
    assert "feature_columns" in schema
    assert schema["feature_columns"]
