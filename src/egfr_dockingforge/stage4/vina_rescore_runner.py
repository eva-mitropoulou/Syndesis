from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

from egfr_dockingforge.common.io import write_table
from egfr_dockingforge.stage4.schemas import EMPIRICAL_SCORE_COLUMNS


def build_empirical_scores(tasks: pd.DataFrame, gnina_scores: pd.DataFrame, paths: dict[str, Path]) -> pd.DataFrame:
    gnina = gnina_scores[["pose_id", "gnina_empirical_affinity"]] if not gnina_scores.empty else pd.DataFrame(columns=["pose_id", "gnina_empirical_affinity"])
    frame = tasks[["pose_id", "original_docking_score"]].merge(gnina, on="pose_id", how="left")
    frame["vina_rescore"] = None
    frame["vinardo_rescore"] = None
    frame["empirical_rescoring_status"] = "original_score_only"
    frame["warnings_json"] = json.dumps(["Vina/Vinardo rescoring not enabled in Stage 4 default config."])
    frame = frame[EMPIRICAL_SCORE_COLUMNS]
    write_table(paths["processed"] / "empirical_scores.parquet", frame)
    write_table(paths["processed"] / "empirical_scores.csv", frame)
    return frame

