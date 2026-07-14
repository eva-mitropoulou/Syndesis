from __future__ import annotations

from pathlib import Path
from typing import Any

import pandas as pd

from syndesis.common.io import ensure_dir, load_yaml, project_root, resolve_path


def load_stage4_config(config_path: str | Path) -> dict[str, Any]:
    path = resolve_path(config_path)
    config = load_yaml(path)
    config["_config_path"] = str(path)
    config["_project_root"] = str(project_root())
    return config


def stage4_paths(config: dict[str, Any]) -> dict[str, Path]:
    root = Path(config["_project_root"])
    paths = {key: resolve_path(value, root) for key, value in config["paths"].items()}
    for path in paths.values():
        ensure_dir(path)
    return paths


def input_path(config: dict[str, Any], key: str) -> Path:
    return resolve_path(config["inputs"][key], Path(config["_project_root"]))


def load_stage3_inputs(config: dict[str, Any]) -> dict[str, pd.DataFrame]:
    return {
        "poses": pd.read_parquet(input_path(config, "docked_poses")),
        "runs": pd.read_parquet(input_path(config, "docking_runs")),
        "rmsd": pd.read_parquet(input_path(config, "pose_rmsd")),
        "sanity": pd.read_parquet(input_path(config, "pose_sanity")),
        "labels": pd.read_parquet(input_path(config, "stage3_pose_labels")),
        "tasks": pd.read_parquet(input_path(config, "docking_task_matrix")),
    }

