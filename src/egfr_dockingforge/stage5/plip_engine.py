from __future__ import annotations

import importlib.util
import json
from pathlib import Path
from typing import Any

import pandas as pd

from egfr_dockingforge.common.io import write_table
from egfr_dockingforge.stage5.schemas import PLIP_CROSSCHECK_COLUMNS


def run_plip_crosscheck(config: dict[str, Any], paths: dict[str, Path]) -> pd.DataFrame:
    enabled = bool(config.get("plip", {}).get("enabled", False))
    installed = importlib.util.find_spec("plip") is not None
    status = "disabled" if not enabled else ("available_not_run_in_mvp" if installed else "plip_not_installed")
    row = {
        "analysis_id": "plip_crosscheck_status",
        "pose_id": None,
        "complex_id": None,
        "entity_type": "stage5",
        "plip_status": status,
        "plip_version": None,
        "plip_interactions_json": json.dumps([]),
        "prolif_interactions_json": json.dumps([]),
        "agreement_score": None,
        "major_disagreements_json": json.dumps([]),
        "disagreement_reason": None,
        "use_for_report_flag": False,
        "warnings_json": json.dumps([] if not enabled else ["plip_crosscheck_is_optional_and_separate"]),
    }
    frame = pd.DataFrame([row], columns=PLIP_CROSSCHECK_COLUMNS)
    write_table(paths["processed"] / "plip_crosscheck.parquet", frame)
    write_table(paths["processed"] / "plip_crosscheck.csv", frame)
    return frame
