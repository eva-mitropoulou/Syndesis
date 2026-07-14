from __future__ import annotations

from pathlib import Path
from typing import Any

import pandas as pd

from syndesis.common.io import ensure_dir, load_yaml, project_root, resolve_path


def load_stage12_config(config_path: str | Path) -> dict[str, Any]:
    config = load_yaml(resolve_path(config_path))
    config["_config_path"] = str(resolve_path(config_path))
    return config


def stage12_paths(config: dict[str, Any]) -> dict[str, Path]:
    outputs = config["outputs"]
    paths = {
        "processed": resolve_path(outputs["processed_dir"]),
        "figures": resolve_path(outputs["figures_dir"]),
        "dossier": resolve_path(outputs["dossier_dir"]),
        "cards": resolve_path(outputs["dossier_dir"]) / "cards",
        "html": resolve_path(outputs["dossier_dir"]) / "html",
        "reports": resolve_path(outputs["reports_dir"]),
        "model_cards": resolve_path(outputs["model_cards_dir"]),
        "dataset_cards": resolve_path(outputs["dataset_cards_dir"]),
        "bundle": resolve_path(outputs["reproducibility_bundle_dir"]),
    }
    for path in paths.values():
        ensure_dir(path)
    return paths


def _read_table(path: str | Path) -> pd.DataFrame | None:
    resolved = resolve_path(path)
    if not resolved.exists():
        return None
    if resolved.suffix == ".parquet":
        return pd.read_parquet(resolved)
    if resolved.suffix == ".csv":
        return pd.read_csv(resolved)
    return None


def load_stage12_inputs(config: dict[str, Any]) -> dict[str, pd.DataFrame | None]:
    return {name: _read_table(path) for name, path in config.get("inputs", {}).items() if name != "stage12_sources"}


def relative_to_root(path: str | Path | None) -> str | None:
    if path in (None, ""):
        return None
    try:
        return str(resolve_path(path).relative_to(project_root()))
    except ValueError:
        return str(path)
