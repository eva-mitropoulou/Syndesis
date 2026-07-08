from __future__ import annotations

from html import escape
from pathlib import Path

import pandas as pd

from egfr_dockingforge.common.io import ensure_dir


def _table(frame: pd.DataFrame, limit: int = 30) -> str:
    if frame.empty:
        return "<p>No rows.</p>"
    return frame.head(limit).to_html(index=False, escape=True)


def write_stage5_report(
    native_long: pd.DataFrame,
    key: pd.DataFrame,
    docked_long: pd.DataFrame,
    recovery: pd.DataFrame,
    clusters: pd.DataFrame,
    final_labels: pd.DataFrame,
    plip: pd.DataFrame,
    out_path: str | Path,
) -> Path:
    out = Path(out_path)
    ensure_dir(out.parent)
    engine = "ProLIF"
    label_counts = final_labels["final_pose_label"].value_counts(dropna=False).rename_axis("final_pose_label").reset_index(name="count") if not final_labels.empty else pd.DataFrame()
    recovery_counts = recovery["interaction_recovery_label"].value_counts(dropna=False).rename_axis("interaction_recovery_label").reset_index(name="count") if not recovery.empty else pd.DataFrame()
    cluster_counts = clusters["cluster_label"].value_counts(dropna=False).rename_axis("cluster_label").reset_index(name="count") if not clusters.empty else pd.DataFrame()
    mismatch_good_rmsd = final_labels[final_labels["final_pose_label"].eq("rmsd_good_interactions_poor")]
    mismatch_bad_rmsd = final_labels[final_labels["final_pose_label"].eq("rmsd_poor_interactions_good")]
    html = f"""<!doctype html>
<html>
<head>
  <meta charset="utf-8">
  <title>Stage 5 Interaction Atlas</title>
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
  <h1>Stage 5 Interaction-Fingerprint Atlas</h1>
  <p class="note">Stage 5 combines native EGFR interaction fingerprints, docked-pose interaction recovery, and RMSD/sanity checks. GNINA scores are retained as features but are not label criteria.</p>
  <h2>Run Summary</h2>
  <ul>
    <li>Native interaction rows: {len(native_long)}</li>
    <li>Docked interaction rows: {len(docked_long)}</li>
    <li>Docked poses with recovery rows: {len(recovery)}</li>
    <li>Key EGFR interactions: {len(key)}</li>
    <li>Interaction engine: {escape(engine)}</li>
  </ul>
  <h2>Key EGFR Interactions</h2>
  {_table(key, 50)}
  <h2>Interaction Recovery Labels</h2>
  {_table(recovery_counts, 50)}
  <h2>Final Pose Labels</h2>
  {_table(label_counts, 50)}
  <h2>Binding-Mode Clusters</h2>
  {_table(cluster_counts, 50)}
  <h2>RMSD Good but Interactions Poor</h2>
  {_table(mismatch_good_rmsd, 30)}
  <h2>RMSD Poor but Interactions Preserved</h2>
  {_table(mismatch_bad_rmsd, 30)}
  <h2>PLIP Cross-Check</h2>
  {_table(plip, 20)}
  <h2>Hydrogen Handling</h2>
  <p>Stage 5 requires ProLIF inputs with explicit hydrogens and fails the run when required hydrogens cannot be prepared.</p>
</body>
</html>
"""
    out.write_text(html, encoding="utf-8")
    return out
