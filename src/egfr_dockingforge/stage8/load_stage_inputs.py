from __future__ import annotations

from pathlib import Path
from typing import Any

import pandas as pd

from egfr_dockingforge.common.io import ensure_dir, load_yaml, project_root, resolve_path


def load_stage8_config(config_path: str | Path) -> dict[str, Any]:
    return load_yaml(resolve_path(config_path))


def stage8_paths(config: dict[str, Any]) -> dict[str, Path]:
    root = project_root()
    processed = resolve_path(config["paths"]["processed"], root)
    paths = {"processed": processed, "reports": resolve_path(config["paths"]["reports"], root)}
    for extra in ["docking_outputs/raw", "docking_outputs/poses", "docking_outputs/logs", "gnina_logs", "prepared_ligands", "prolif_ligands", "prolif_receptors"]:
        paths[extra.replace("/", "_")] = ensure_dir(processed / extra)
    for path in paths.values():
        ensure_dir(path)
    return paths


def read_table(path: str | Path) -> pd.DataFrame:
    resolved = resolve_path(path)
    if not resolved.exists():
        raise FileNotFoundError(f"Missing Stage 8 input: {resolved}")
    return pd.read_csv(resolved) if resolved.suffix == ".csv" else pd.read_parquet(resolved)


def load_stage8_inputs(config: dict[str, Any]) -> dict[str, pd.DataFrame]:
    return {name: read_table(path) for name, path in config["inputs"].items() if str(path).endswith((".parquet", ".csv"))}
