from __future__ import annotations

from pathlib import Path
from typing import Any

import pandas as pd

from egfr_dockingforge.common.io import ensure_dir, load_yaml, project_root, resolve_path


def load_stage9_config(config_path: str | Path) -> dict[str, Any]:
    return load_yaml(resolve_path(config_path))


def stage9_paths(config: dict[str, Any]) -> dict[str, Path]:
    root = project_root()
    processed = ensure_dir(resolve_path(config["paths"]["processed"], root))
    reports = ensure_dir(resolve_path(config["paths"]["reports"], root))
    analog_sdf = ensure_dir(resolve_path(config["paths"]["analog_sdf_dir"], root))
    return {"processed": processed, "reports": reports, "analog_sdf_dir": analog_sdf}


def read_table(path: str | Path) -> pd.DataFrame:
    resolved = resolve_path(path)
    if not resolved.exists():
        raise FileNotFoundError(f"Missing Stage 9 input: {resolved}")
    return pd.read_csv(resolved) if resolved.suffix == ".csv" else pd.read_parquet(resolved)


def load_stage9_inputs(config: dict[str, Any]) -> dict[str, pd.DataFrame]:
    return {
        name: read_table(path)
        for name, path in config["inputs"].items()
        if str(path).endswith((".parquet", ".csv"))
    }
