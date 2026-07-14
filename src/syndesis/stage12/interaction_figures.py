from __future__ import annotations

from pathlib import Path

import pandas as pd

from syndesis.stage12.structure_figures import _svg_text


def write_interaction_recovery_distribution(selection: pd.DataFrame, out_dir: Path) -> Path:
    values = selection["best_key_interaction_recall_consensus"].fillna(0).round(3).tolist()
    return _svg_text(out_dir / "prolif_interaction_recovery_distribution.svg", "ProLIF interaction recovery distribution", [f"recall values: {values}"])


def write_evidence_matrix(selection: pd.DataFrame, out_dir: Path) -> Path:
    lines = [
        f"{r.final_candidate_id}: score={r.final_candidate_score:.3f}, pose={r.best_pose_confidence}, md={r.md_stability_label_if_available}, decision={r.decision_label}"
        for r in selection.itertuples()
    ]
    return _svg_text(out_dir / "final_candidate_evidence_matrix.svg", "Final candidate evidence matrix", lines, width=1080, height=620)
