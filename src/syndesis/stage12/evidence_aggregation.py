from __future__ import annotations

from pathlib import Path
from typing import Any

import pandas as pd

from syndesis.stage12.load_stage_inputs import load_stage12_config, load_stage12_inputs, stage12_paths


def load_selection_and_inputs(config_path: str | Path) -> tuple[dict[str, Any], dict[str, Path], dict[str, pd.DataFrame | None], pd.DataFrame]:
    config = load_stage12_config(config_path)
    paths = stage12_paths(config)
    selection_path = paths["processed"] / "final_candidate_selection.parquet"
    if not selection_path.exists():
        from syndesis.stage12.final_candidate_selection import build_final_candidate_table

        build_final_candidate_table(config_path)
    selection = pd.read_parquet(selection_path)
    return config, paths, load_stage12_inputs(config), selection


def best_pose_file(pose_id: str | None) -> str | None:
    if not pose_id:
        return None
    candidates = [
        Path("data/processed/stage8/docking_outputs/poses") / f"{pose_id}.pdbqt",
        Path("data/processed/stage8/gnina_ligands") / f"{pose_id}.pdb",
        Path("data/processed/stage8/prolif_ligands") / f"{pose_id}.pose_template.sdf",
    ]
    for path in candidates:
        if path.exists():
            return str(path)
    return None
