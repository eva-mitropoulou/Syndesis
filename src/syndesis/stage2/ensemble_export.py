from __future__ import annotations

import json
import shutil
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from syndesis.common.io import ensure_dir
from syndesis.stage2.pocket_features import pdb_heavy_atom_coords
from syndesis.stage2.schemas import ENSEMBLE_COLUMNS


def docking_box(coords: np.ndarray, config: dict[str, Any]) -> tuple[list[float], list[float]]:
    box_cfg = config["docking_box_suggestions"]
    if coords.size == 0:
        return [None, None, None], [None, None, None]
    center = coords.mean(axis=0)
    span = coords.max(axis=0) - coords.min(axis=0)
    size = span + float(box_cfg.get("padding_angstrom", 8.0))
    size = np.maximum(size, float(box_cfg.get("min_side_angstrom", 18.0)))
    size = np.minimum(size, float(box_cfg.get("max_side_angstrom", 30.0)))
    return [round(float(x), 3) for x in center], [round(float(x), 3) for x in size]


def export_ensemble(selected: pd.DataFrame, aligned_lookup: pd.DataFrame, out_dir: Path, config: dict[str, Any]) -> pd.DataFrame:
    ensure_dir(out_dir)
    for stale in out_dir.glob("*.pdb"):
        stale.unlink()
    aligned = aligned_lookup.set_index("receptor_id")["aligned_receptor_file_path"].to_dict()
    rows = []
    for _, row in selected.iterrows():
        rid = row["receptor_id"]
        target = out_dir / f"{rid}.pdb"
        shutil.copyfile(row["receptor_file_path"], target)
        coords = pdb_heavy_atom_coords(row["native_ligand_sdf_path"])
        center, size = docking_box(coords, config)
        native_centroid = [
            row["native_ligand_centroid_x"],
            row["native_ligand_centroid_y"],
            row["native_ligand_centroid_z"],
        ]
        record = row.to_dict()
        record.update(
            {
                "receptor_file_path": str(target),
                "aligned_receptor_file_path": aligned.get(rid),
                "native_ligand_centroid": json.dumps(native_centroid),
                "suggested_docking_box_center": json.dumps(center),
                "suggested_docking_box_size": json.dumps(size),
                "stage3_validation_status": "pending",
            }
        )
        rows.append({column: record.get(column) for column in ENSEMBLE_COLUMNS})
    return pd.DataFrame(rows, columns=ENSEMBLE_COLUMNS)
