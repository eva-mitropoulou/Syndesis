from __future__ import annotations

from pathlib import Path

import pandas as pd

from syndesis.common.io import ensure_dir


def counts(frame: pd.DataFrame, column: str) -> str:
    if frame.empty or column not in frame:
        return "<p>No data.</p>"
    return frame[column].astype("string").fillna("unknown").value_counts().reset_index().to_html(index=False)


def table(frame: pd.DataFrame, columns: list[str], limit: int = 50) -> str:
    if frame.empty:
        return "<p>No rows.</p>"
    return frame[[c for c in columns if c in frame.columns]].head(limit).to_html(index=False)


def write_stage3_report(tasks: pd.DataFrame, runs: pd.DataFrame, metrics: pd.DataFrame, validation: pd.DataFrame, out_path: Path) -> Path:
    ensure_dir(out_path.parent)
    html = f"""<!doctype html>
<html><head><meta charset="utf-8"><title>Syndesis Stage 3</title></head>
<body>
<h1>Stage 3 Redocking And Cross-docking Benchmark</h1>
<h2>Summary</h2>
<ul>
<li>Total docking tasks: {len(tasks)}</li>
<li>Successful runs: {(runs['status'] == 'success').sum() if not runs.empty else 0}</li>
<li>Failed/skipped runs: {(runs['status'] != 'success').sum() if not runs.empty else 0}</li>
</ul>
<h2>Task Types</h2>{counts(tasks, "task_type")}
<h2>Run Status</h2>{counts(runs, "status")}
<h2>Failure Categories</h2>{counts(metrics, "failure_category")}
<h2>Receptor Recommendations</h2>{table(validation, ["receptor_id", "pdb_id", "dominant_failure_mode", "keep_for_stage4_flag", "recommendation_reason"])}
<h2>Controls</h2>{table(tasks[tasks["native_receptor_id"].str.contains("1m17|1xkk", case=False, regex=True)], ["docking_task_id", "task_type", "native_receptor_id", "target_receptor_id"], 20)}
<h2>Scope</h2><p>Stage 3 reports docking geometry, RMSD, and receptor-level pose reproducibility. Later-stage evidence is reported in the corresponding stage reports.</p>
</body></html>
"""
    out_path.write_text(html, encoding="utf-8")
    return out_path
