"""Audit EF1% cutoff ties for the frozen paper-ranking scores."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from submission_robustness_analysis import TARGETS, load_target  # noqa: E402
from syndesis.enrichment.metrics import enrichment_factor  # noqa: E402


def audit_arm(labels: np.ndarray, scores: np.ndarray, *, seed: int = 807) -> dict[str, float | int]:
    cutoff = max(1, round(len(labels) * 0.01))
    stable_order = np.argsort(-scores, kind="mergesort")
    boundary = scores[stable_order[cutoff - 1]]
    above = scores > boundary
    tied = scores == boundary
    slots = cutoff - int(above.sum())
    tied_labels = labels[tied]
    stable_ef = enrichment_factor(labels, scores, 0.01)
    fractional_top_actives = labels[above].sum() + slots * tied_labels.mean()
    fractional_ef = fractional_top_actives / cutoff / labels.mean()

    random_efs = []
    rng = np.random.default_rng(seed)
    for _ in range(1000):
        tie_indices = np.flatnonzero(tied)
        chosen = rng.choice(tie_indices, size=slots, replace=False)
        active_count = int(labels[above].sum() + labels[chosen].sum())
        random_efs.append(active_count / cutoff / labels.mean())

    return {
        "n_ligands": len(labels),
        "top_fraction": 0.01,
        "cutoff_n": cutoff,
        "boundary_score": boundary,
        "n_strictly_above_boundary": int(above.sum()),
        "boundary_tie_size": int(tied.sum()),
        "boundary_tie_active_count": int(tied_labels.sum()),
        "stable_order_ef1": stable_ef,
        "fractional_tie_ef1": fractional_ef,
        "random_tie_ef1_min": float(np.min(random_efs)),
        "random_tie_ef1_max": float(np.max(random_efs)),
    }


def main() -> None:
    rows = []
    for target, config in TARGETS.items():
        _, per_ligand, _ = load_target(config)
        labels = per_ligand["label"].to_numpy(int)
        for arm in ("gnina", "coupled"):
            rows.append({"target": target, "arm": arm, **audit_arm(labels, per_ligand[arm].to_numpy(float))})
    out = ROOT / "results" / "robustness" / "ef1_tie_handling_audit.csv"
    pd.DataFrame(rows).to_csv(out, index=False)
    print(out)


if __name__ == "__main__":
    main()
