from __future__ import annotations

from pathlib import Path

import pandas as pd


def _table(df: pd.DataFrame, n: int = 20) -> str:
    return df.head(n).to_html(index=False, escape=True)


def write_stage8_report(paths: dict) -> Path:
    report = paths["reports"] / "08_candidate_screening.html"
    manifest = pd.read_parquet(paths["processed"] / "screening_candidate_manifest.parquet")
    tasks = pd.read_parquet(paths["processed"] / "screening_task_matrix.parquet")
    runs = pd.read_parquet(paths["processed"] / "screening_docking_runs.parquet")
    ranked = pd.read_parquet(paths["processed"] / "ranked_candidates.parquet")
    diag = pd.read_parquet(paths["processed"] / "known_control_diagnostics.parquet")
    html = f"""<!doctype html><html><head><meta charset="utf-8"><title>Stage 8 Candidate Screening</title>
<style>body{{font-family:Arial,sans-serif;margin:32px}}table{{border-collapse:collapse}}td,th{{border:1px solid #ccc;padding:4px 8px}}</style></head><body>
<h1>Stage 8 Source-aware Structure-based Candidate Screening</h1>
<p>Screening manifest rows: {len(manifest)}. Tasks: {len(tasks)}. Docking runs: {len(runs)}. Ranked molecules: {len(ranked)}.</p>
<h2>Source Composition</h2>{_table(manifest['source'].value_counts().rename_axis('source').reset_index(name='count'))}
<h2>Docking Completion</h2>{_table(runs['status'].value_counts().rename_axis('status').reset_index(name='count'))}
<h2>Known-control Diagnostics</h2>{_table(diag)}
<h2>Top Ranked Candidates</h2>{_table(ranked, 50)}
<h2>Non-claims</h2><p>Known-control diagnostics are reported separately from prospective candidate ranking. No experimental activity is claimed.</p>
</body></html>"""
    report.write_text(html, encoding="utf-8")
    return report
