from __future__ import annotations

from html import escape
from pathlib import Path

import pandas as pd

from egfr_dockingforge.common.io import ensure_dir


def value_counts_table(frame: pd.DataFrame, column: str) -> str:
    if frame.empty or column not in frame:
        return "<p>No data.</p>"
    counts = frame[column].astype("string").fillna("unknown").value_counts().reset_index()
    counts.columns = [column, "count"]
    return counts.to_html(index=False, escape=True)


def small_table(frame: pd.DataFrame, columns: list[str], limit: int = 25) -> str:
    if frame.empty:
        return "<p>No rows.</p>"
    available = [column for column in columns if column in frame.columns]
    return frame[available].head(limit).to_html(index=False, escape=True)


def write_stage1_report(benchmark: pd.DataFrame, rejected: pd.DataFrame, out_path: Path) -> Path:
    ensure_dir(out_path.parent)
    included = benchmark[benchmark["include_in_stage1_benchmark"] == True] if not benchmark.empty else benchmark
    html = f"""<!doctype html>
<html>
<head><meta charset="utf-8"><title>EGFR DockingForge Stage 1</title></head>
<body>
<h1>EGFR Co-crystal Benchmark</h1>
<h2>Summary</h2>
<ul>
<li>Candidate chain-ligand complexes: {len(benchmark)}</li>
<li>Retained complexes: {len(included)}</li>
<li>Rejected complexes: {len(rejected)}</li>
</ul>
<h2>Rejection Reasons</h2>
{value_counts_table(rejected, "exclusion_reason")}
<h2>Quality Tier Distribution</h2>
{value_counts_table(benchmark, "quality_tier")}
<h2>Resolution Distribution</h2>
{benchmark["resolution_angstrom"].describe().to_frame().to_html(escape=True) if "resolution_angstrom" in benchmark else "<p>No data.</p>"}
<h2>Ligand Heavy Atom Distribution</h2>
{benchmark["ligand_heavy_atom_count"].describe().to_frame().to_html(escape=True) if "ligand_heavy_atom_count" in benchmark else "<p>No data.</p>"}
<h2>Receptor-state Distribution</h2>
{value_counts_table(benchmark, "kincore_state")}
<h2>Mutation-status Distribution</h2>
{value_counts_table(benchmark, "mutation_flag")}
<h2>Active-site Completeness</h2>
{benchmark["active_site_completeness_score"].describe().to_frame().to_html(escape=True) if "active_site_completeness_score" in benchmark else "<p>No data.</p>"}
<h2>Tier A Complexes</h2>
{small_table(benchmark[benchmark["quality_tier"] == "Tier A"], ["complex_id", "pdb_id", "ligand_comp_id", "resolution_angstrom", "quality_score"])}
<h2>Tier B Complexes</h2>
{small_table(benchmark[benchmark["quality_tier"] == "Tier B"], ["complex_id", "pdb_id", "ligand_comp_id", "resolution_angstrom", "quality_score"])}
<h2>Rejected Complexes</h2>
{small_table(rejected, ["complex_id", "pdb_id", "ligand_comp_id", "exclusion_reason"], limit=100)}
<h2>Source and Provenance Summary</h2>
<p>Raw RCSB mmCIF, metadata, validation, CCD, and optional KLIFS files are stored under <code>data/raw/stage1</code>.</p>
<h2>Warning Summary</h2>
{small_table(benchmark[benchmark["warnings_json"].fillna("[]") != "[]"], ["complex_id", "warnings_json"], limit=100)}
</body>
</html>
"""
    out_path.write_text(html, encoding="utf-8")
    return out_path
