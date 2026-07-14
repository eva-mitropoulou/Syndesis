from __future__ import annotations

from pathlib import Path
from typing import Any

import pandas as pd

from syndesis.common.io import ensure_dir, load_yaml, project_root, resolve_path


def load_stage7_config(config_path: str | Path) -> dict[str, Any]:
    return load_yaml(resolve_path(config_path))


def stage7_paths(config: dict[str, Any]) -> dict[str, Path]:
    root = project_root()
    paths = {
        "raw": resolve_path(config["paths"]["raw"], root),
        "processed": resolve_path(config["paths"]["processed"], root),
        "reports": resolve_path(config["paths"]["reports"], root),
    }
    for path in paths.values():
        ensure_dir(path)
    ensure_dir(paths["processed"] / "prepared_ligands")
    ensure_dir(paths["raw"] / "chembl")
    ensure_dir(paths["raw"] / "bindingdb")
    ensure_dir(paths["raw"] / "pubchem")
    ensure_dir(paths["raw"] / "vendor")
    return paths


def read_table(path: str | Path) -> pd.DataFrame:
    resolved = resolve_path(path)
    if not resolved.exists():
        raise FileNotFoundError(f"Missing Stage 7 input: {resolved}")
    if resolved.suffix == ".csv":
        return pd.read_csv(resolved)
    return pd.read_parquet(resolved)


def load_stage7_inputs(config: dict[str, Any]) -> dict[str, pd.DataFrame]:
    return {name: read_table(path) for name, path in config["inputs"].items() if not str(path).endswith(".json")}
