from __future__ import annotations

from pathlib import Path
from typing import Any

import pandas as pd

from egfr_dockingforge.common.io import ensure_dir, load_yaml, project_root, resolve_path


def load_stage11_config(config_path: str | Path) -> dict[str, Any]:
    return load_yaml(resolve_path(config_path))


def stage11_paths(config: dict[str, Any]) -> dict[str, Path]:
    root = project_root()
    return {
        "processed": ensure_dir(resolve_path(config["paths"]["processed"], root)),
        "reports": ensure_dir(resolve_path(config["paths"]["reports"], root)),
        "md_root": ensure_dir(resolve_path(config["paths"]["md_root"], root)),
        "user_cgenff_str_dir": ensure_dir(resolve_path(config["paths"]["user_cgenff_str_dir"], root)),
    }


def read_optional(path: str | Path) -> pd.DataFrame | None:
    p = resolve_path(path)
    if not p.exists():
        return None
    return pd.read_csv(p) if p.suffix == ".csv" else pd.read_parquet(p)


def load_stage11_inputs(config: dict[str, Any]) -> dict[str, pd.DataFrame | None]:
    return {name: read_optional(path) for name, path in config["inputs"].items()}
