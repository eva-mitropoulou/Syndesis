from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

from egfr_dockingforge.common.io import write_table
from egfr_dockingforge.stage3.schemas import RECEPTOR_VALIDATION_COLUMNS


def receptor_validation_summary(ensemble: pd.DataFrame, metrics: pd.DataFrame, paths: dict[str, Path]) -> pd.DataFrame:
    rows = []
    for _, receptor in ensemble.iterrows():
        rid = receptor["receptor_id"]
        target_metrics = metrics[metrics["target_receptor_id"] == rid] if not metrics.empty else pd.DataFrame()
        redock = target_metrics[target_metrics["task_type"] == "redocking_native_receptor"] if not target_metrics.empty else pd.DataFrame()
        cross = target_metrics[target_metrics["task_type"] != "redocking_native_receptor"] if not target_metrics.empty else pd.DataFrame()
        same = cross[cross["state_match_flag"] == True] if not cross.empty else pd.DataFrame()
        other = cross[cross["state_match_flag"] != True] if not cross.empty else pd.DataFrame()
        def rate(frame: pd.DataFrame, column: str) -> float | None:
            return None if frame.empty else float((frame[column] == True).mean())
        if target_metrics.empty:
            dominant = "no_tasks"
        else:
            dominant = target_metrics["failure_category"].mode().iloc[0]
        keep = dominant not in {"engine_unavailable", "no_tasks"}
        rows.append({
            "receptor_id": rid,
            "pdb_id": receptor["pdb_id"],
            "receptor_state": receptor["state_stratum"],
            "num_redocking_tasks": int(len(redock)),
            "num_crossdocking_tasks": int(len(cross)),
            "redocking_top1_success_rate_strict": rate(redock, "ranking_success_strict_top1"),
            "redocking_top20_success_rate_strict": rate(redock, "sampling_success_strict_top20"),
            "same_state_crossdock_top1_success_rate_strict": rate(same, "ranking_success_strict_top1"),
            "same_state_crossdock_top20_success_rate_strict": rate(same, "sampling_success_strict_top20"),
            "other_state_crossdock_top20_success_rate_relaxed": rate(other, "sampling_success_relaxed_top20"),
            "physical_sanity_pass_rate": None if target_metrics.empty else float(target_metrics["physical_sanity_pass_rate"].dropna().mean()),
            "dominant_failure_mode": dominant,
            "keep_for_stage4_flag": keep,
            "prune_recommendation_flag": False if keep else True,
            "recommendation_reason": "Review receptor performance before pruning; Stage 3 gives recommendations only.",
            "warnings_json": json.dumps(["Do not prune automatically; Stage 3 validation incomplete."]),
        })
    frame = pd.DataFrame(rows, columns=RECEPTOR_VALIDATION_COLUMNS)
    write_table(paths["processed"] / "receptor_stage3_validation.parquet", frame)
    write_table(paths["processed"] / "receptor_stage3_validation.csv", frame)
    return frame
