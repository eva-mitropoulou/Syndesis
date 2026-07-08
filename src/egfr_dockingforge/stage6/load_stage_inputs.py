from __future__ import annotations

from pathlib import Path
from typing import Any

import pandas as pd

from egfr_dockingforge.common.io import ensure_dir, load_yaml, project_root, resolve_path


def load_stage6_config(config_path: str | Path) -> dict[str, Any]:
    return load_yaml(resolve_path(config_path))


def stage6_paths(config: dict[str, Any]) -> dict[str, Path]:
    root = project_root()
    paths = {
        "processed": resolve_path(config["paths"]["processed"], root),
        "models": resolve_path(config["paths"]["models"], root),
        "reports": resolve_path(config["paths"]["reports"], root),
    }
    for path in paths.values():
        ensure_dir(path)
    return paths


def read_table(path: str | Path) -> pd.DataFrame:
    resolved = resolve_path(path)
    if not resolved.exists():
        raise FileNotFoundError(f"Required Stage 6 input is missing: {resolved}")
    if resolved.suffix == ".csv":
        return pd.read_csv(resolved)
    return pd.read_parquet(resolved)


def load_stage6_inputs(config: dict[str, Any]) -> dict[str, pd.DataFrame]:
    return {name: read_table(path) for name, path in config["inputs"].items()}
