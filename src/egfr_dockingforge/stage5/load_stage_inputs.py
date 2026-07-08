from __future__ import annotations

from pathlib import Path
from typing import Any

import pandas as pd

from egfr_dockingforge.common.io import ensure_dir, load_yaml, project_root, resolve_path


def load_stage5_config(config_path: str | Path) -> dict[str, Any]:
    return load_yaml(resolve_path(config_path))


def stage5_paths(config: dict[str, Any]) -> dict[str, Path]:
    root = project_root()
    paths = {
        "processed": resolve_path(config["paths"]["processed"], root),
        "reports": resolve_path(config["paths"]["reports"], root),
    }
    for path in paths.values():
        ensure_dir(path)
    return paths


def _read_table(path: str | Path) -> pd.DataFrame:
    resolved = resolve_path(path)
    if resolved.suffix == ".csv":
        return pd.read_csv(resolved)
    return pd.read_parquet(resolved)


def load_stage5_inputs(config: dict[str, Any]) -> dict[str, pd.DataFrame]:
    inputs = config["inputs"]
    return {
        "benchmark": _read_table(inputs["cocrystal_benchmark"]),
        "ensemble": _read_table(inputs["receptor_ensemble"]),
        "pocket_mapping": _read_table(inputs["pocket_residue_mapping"]),
        "docked_poses": _read_table(inputs["docked_poses"]),
        "labels": _read_table(inputs["stage3_pose_labels"]),
        "pose_scores": _read_table(inputs["pose_score_table"]),
        "stage4_features": _read_table(inputs["stage4_features"]),
    }
