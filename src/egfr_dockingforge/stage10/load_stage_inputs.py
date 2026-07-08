from __future__ import annotations

from pathlib import Path
from typing import Any

import pandas as pd

from egfr_dockingforge.common.io import ensure_dir, load_yaml, project_root, resolve_path


def load_stage10_config(config_path: str | Path) -> dict[str, Any]:
    return load_yaml(resolve_path(config_path))


def stage10_paths(config: dict[str, Any]) -> dict[str, Path]:
    root = project_root()
    return {
        "processed": ensure_dir(resolve_path(config["paths"]["processed"], root)),
        "reports": ensure_dir(resolve_path(config["paths"]["reports"], root)),
        "figures": ensure_dir(resolve_path(config["paths"]["figures"], root)),
    }


def read_optional_table(path: str | Path) -> pd.DataFrame | None:
    resolved = resolve_path(path)
    if not resolved.exists():
        return None
    return pd.read_csv(resolved) if resolved.suffix == ".csv" else pd.read_parquet(resolved)


def load_stage10_inputs(config: dict[str, Any]) -> dict[str, pd.DataFrame | None]:
    return {name: read_optional_table(path) for name, path in config["inputs"].items()}
