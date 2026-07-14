from __future__ import annotations

from html import escape
from pathlib import Path
from typing import Any

import pandas as pd

from syndesis.common.io import ensure_dir


def _table(frame: pd.DataFrame, limit: int = 20) -> str:
    if frame.empty:
        return "<p>No rows.</p>"
    return frame.head(limit).to_html(index=False, escape=True)


def write_stage4_report(
    tasks: pd.DataFrame,
    gnina: pd.DataFrame,
    metrics: pd.DataFrame,
    summary: pd.DataFrame,
    failures: pd.DataFrame,
    raw_runs: pd.DataFrame,
    out_path: str | Path,
) -> Path:
    out = Path(out_path)
    ensure_dir(out.parent)
    version = ""
    model = ""
    if not gnina.empty:
        version = str(gnina["gnina_version"].dropna().iloc[0]) if gnina["gnina_version"].notna().any() else ""
        model = str(gnina["gnina_model"].dropna().iloc[0]) if gnina["gnina_model"].notna().any() else ""
    scored = int(gnina["rescoring_status"].eq("success").sum()) if not gnina.empty else 0
    failed = int(gnina["rescoring_status"].ne("success").sum()) if not gnina.empty else 0
    runtime_total = float(pd.to_numeric(raw_runs.get("runtime_seconds"), errors="coerce").sum()) if not raw_runs.empty else 0.0
    task_counts = tasks["task_status"].value_counts(dropna=False).rename_axis("task_status").reset_index(name="count")
    label_counts = gnina["stage3_pose_label"].value_counts(dropna=False).rename_axis("stage3_pose_label").reset_index(name="count") if not gnina.empty else pd.DataFrame()
    html = f"""<!doctype html>
<html>
<head>
  <meta charset="utf-8">
  <title>Stage 4 ML Rescoring</title>
  <style>
    body {{ font-family: Arial, sans-serif; margin: 32px; color: #202124; }}
    h1, h2 {{ color: #1f2933; }}
    table {{ border-collapse: collapse; margin: 12px 0 24px 0; font-size: 13px; }}
    th, td {{ border: 1px solid #d6d9de; padding: 6px 8px; vertical-align: top; }}
    th {{ background: #eef2f6; }}
    .note {{ background: #f8fafc; border-left: 4px solid #6b7280; padding: 10px 12px; }}
  </style>
</head>
<body>
  <h1>Stage 4 ML-Based Pose Rescoring</h1>
  <p class="note">GNINA scores are treated as pose-confidence features and are not used alone for candidate acceptance.</p>
  <h2>Run Summary</h2>
  <ul>
    <li>Rescoring tasks: {len(tasks)}</li>
    <li>GNINA scored poses: {scored}</li>
    <li>GNINA failed/missing poses: {failed}</li>
    <li>GNINA version: {escape(version)}</li>
    <li>GNINA model: {escape(model)}</li>
    <li>Total GNINA runtime seconds: {runtime_total:.1f}</li>
  </ul>
  <h2>Task Status</h2>
  {_table(task_counts)}
  <h2>Stage 3 Label Distribution Among Scored Poses</h2>
  {_table(label_counts)}
  <h2>Scorer Comparison Summary</h2>
  {_table(summary, 50)}
  <h2>Per-Task Metrics</h2>
  {_table(metrics, 50)}
  <h2>Failure Examples</h2>
  {_table(failures, 50)}
  <h2>Recommendation</h2>
  <p>Use GNINA CNNscore, CNNaffinity, and empirical affinity as Stage 6 features alongside original docking scores. Do not use GNINA alone for candidate acceptance.</p>
</body>
</html>
"""
    out.write_text(html, encoding="utf-8")
    return out
