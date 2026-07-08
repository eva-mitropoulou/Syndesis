from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

from egfr_dockingforge.common.io import write_table


def compute_diversity_novelty(master: pd.DataFrame, paths: dict[str, Path]) -> pd.DataFrame:
    rows = []
    for (sid, seed), g in master.groupby(["strategy_id", "seed_id"], dropna=False):
        unique_scaffolds = g["scaffold_id"].nunique()
        rows.append(
            {
                "strategy_id": sid,
                "seed_id": seed,
                "internal_diversity": float(1 - g["parent_tanimoto"].mean()) if len(g) else 0.0,
                "unique_scaffold_count": int(unique_scaffolds),
                "unique_scaffold_rate": unique_scaffolds / max(len(g), 1),
                "mean_parent_tanimoto": float(g["parent_tanimoto"].mean()) if len(g) else 0.0,
                "median_parent_tanimoto": float(g["parent_tanimoto"].median()) if len(g) else 0.0,
                "known_duplicate_rate": float((g["unique_flag"] == False).mean()) if len(g) else 0.0,
                "close_analog_rate": float((g["parent_tanimoto"] >= 0.7).mean()) if len(g) else 0.0,
                "scaffold_novel_rate": 0.0,
                "mode_collapse_flag": bool(unique_scaffolds <= 1 and len(g) > 5),
                "warnings_json": json.dumps([]),
            }
        )
    out = pd.DataFrame(rows)
    write_table(paths["processed"] / "diversity_novelty_metrics.parquet", out)
    write_table(paths["processed"] / "diversity_novelty_metrics.csv", out)
    return out
