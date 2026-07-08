from __future__ import annotations

from pathlib import Path

import pandas as pd

from egfr_dockingforge.common.io import ensure_dir


def counts(frame: pd.DataFrame, column: str) -> str:
    if frame.empty or column not in frame:
        return "<p>No data.</p>"
    return frame[column].astype("string").fillna("unknown").value_counts().reset_index().to_html(index=False)


def table(frame: pd.DataFrame, columns: list[str], limit: int = 50) -> str:
    if frame.empty:
        return "<p>No rows.</p>"
    return frame[[c for c in columns if c in frame.columns]].head(limit).to_html(index=False)


def write_stage2_report(
    features: pd.DataFrame,
    excluded: pd.DataFrame,
    clusters: pd.DataFrame,
    ensemble: pd.DataFrame,
    holdout: pd.DataFrame,
    out_path: Path,
) -> Path:
    ensure_dir(out_path.parent)
    html = f"""<!doctype html>
<html><head><meta charset="utf-8"><title>EGFR DockingForge Stage 2</title></head>
<body>
<h1>EGFR Receptor Ensemble v1</h1>
<h2>Summary</h2>
<ul>
<li>Stage 2 passing receptors: {len(features)}</li>
<li>Excluded receptors: {len(excluded)}</li>
<li>Selected ensemble receptors: {len(ensemble)}</li>
<li>Holdout receptors: {len(holdout)}</li>
</ul>
<h2>State Distribution</h2>{counts(features, "state_stratum")}
<h2>Mutation Distribution</h2>{counts(features, "mutation_flag")}
<h2>Quality Tier Distribution</h2>{counts(features, "quality_tier")}
<h2>Pocket Completeness</h2>{features["active_site_completeness_score"].describe().to_frame().to_html() if not features.empty else "<p>No data.</p>"}
<h2>KLIFS/KinCore Metadata Coverage</h2>{counts(features.assign(kincore_present=features["kincore_activity_label"].notna()), "kincore_present") if not features.empty else "<p>No data.</p>"}
<h2>Cluster Table</h2>{table(clusters, ["receptor_id", "state_stratum", "cluster_id", "cluster_medoid_flag", "nearest_neighbor_distance"])}
<h2>Selected Ensemble</h2>{table(ensemble, ["receptor_id", "pdb_id", "ligand_comp_id", "selected_role", "selected_reason", "state_stratum", "cluster_id", "suggested_docking_box_center", "suggested_docking_box_size"])}
<h2>Holdout Receptors</h2>{table(holdout, ["receptor_id", "pdb_id", "ligand_comp_id", "state_stratum", "cluster_id", "quality_tier"])}
<h2>Excluded Receptors</h2>{table(excluded, ["complex_id", "pdb_id", "ligand_comp_id", "stage2_exclusion_reason"], limit=100)}
<h2>Warnings</h2>{table(features[features["warnings_json"].fillna("[]") != "[]"], ["receptor_id", "warnings_json"], limit=100) if not features.empty else "<p>No data.</p>"}
</body></html>
"""
    out_path.write_text(html, encoding="utf-8")
    return out_path

