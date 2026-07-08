from __future__ import annotations

import importlib.util
import json
from pathlib import Path
from typing import Any

import pandas as pd

from egfr_dockingforge.common.io import write_table


def run_reinvent4_baseline(config: dict[str, Any], paths: dict[str, Path]) -> pd.DataFrame:
    installed = importlib.util.find_spec("reinvent") is not None
    row = {
        "strategy_name": "reinvent4_baseline",
        "status": "not_run_reinvent4_not_installed" if not installed else "not_run_requires_reinvent4_runfile",
        "num_raw_proposals": 0,
        "warnings_json": json.dumps([] if installed else ["reinvent4_python_package_not_found"]),
    }
    out = pd.DataFrame([row])
    write_table(paths["processed"] / "reinvent4_baseline_status.parquet", out)
    write_table(paths["processed"] / "reinvent4_baseline_status.csv", out)
    return out
