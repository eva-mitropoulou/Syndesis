from __future__ import annotations

from pathlib import Path

import pandas as pd

from syndesis.common.io import write_table


def compute_cost_metrics(strategy_metrics: pd.DataFrame, budget: pd.DataFrame, paths: dict[str, Path]) -> pd.DataFrame:
    rows = []
    for row in strategy_metrics.to_dict("records"):
        b = budget[budget["strategy_id"].eq(row["strategy_id"])]
        accepted = max(int(row["num_pre_md_accepted"]), 1)
        gpu = float(b["gpu_seconds"].sum()) if not b.empty else 0.0
        wall = float(b["walltime_seconds"].sum()) if not b.empty else 0.0
        rows.append(
            {
                "strategy_id": row["strategy_id"],
                "cpu_seconds_total": float(b["cpu_seconds"].sum()) if not b.empty else 0.0,
                "gpu_seconds_total": gpu,
                "walltime_seconds_total": wall,
                "num_docking_tasks": int(b["num_docking_tasks"].sum()) if not b.empty else 0,
                "num_gnina_tasks": int(b["num_gnina_tasks"].sum()) if not b.empty else 0,
                "num_prolif_tasks": int(b["num_prolif_tasks"].sum()) if not b.empty else 0,
                "num_md_tasks_if_available": int(b["num_md_tasks_if_available"].sum()) if not b.empty else 0,
                "accepted_analogs_per_gpu_hour": row["num_pre_md_accepted"] / (gpu / 3600) if gpu else 0.0,
                "accepted_analogs_per_wall_hour": row["num_pre_md_accepted"] / (wall / 3600) if wall else 0.0,
                "notes": "cost table derived from Stage 9 logs; no new compute launched in Stage 10",
            }
        )
    out = pd.DataFrame(rows)
    write_table(paths["processed"] / "compute_cost_metrics.parquet", out)
    write_table(paths["processed"] / "compute_cost_metrics.csv", out)
    return out
