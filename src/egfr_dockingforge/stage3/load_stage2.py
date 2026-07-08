from __future__ import annotations

from pathlib import Path
from typing import Any

import pandas as pd

from egfr_dockingforge.common.io import load_yaml, project_root, resolve_path


def load_stage3_config(config_path: str | Path) -> dict[str, Any]:
    return load_yaml(resolve_path(config_path, project_root()))


def stage3_paths(config: dict[str, Any]) -> dict[str, Path]:
    root = project_root()
    return {key: resolve_path(value, root) for key, value in config["paths"].items()}


def load_stage2_ensemble(config: dict[str, Any]) -> pd.DataFrame:
    return pd.read_parquet(resolve_path(config["inputs"]["receptor_ensemble"], project_root()))


def load_stage2_holdout(config: dict[str, Any]) -> pd.DataFrame:
    path = resolve_path(config["inputs"]["receptor_holdout"], project_root())
    return pd.read_parquet(path) if path.exists() else pd.DataFrame()

