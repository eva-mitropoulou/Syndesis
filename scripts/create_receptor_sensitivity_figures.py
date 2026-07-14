#!/usr/bin/env python3
"""Render manuscript receptor-exclusion figures from fixed reported results."""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import Normalize


ROOT = Path(__file__).resolve().parents[1]
OUTS = (ROOT / "figures" / "manuscript", ROOT / "manuscript" / "figures")


def leave_one_receptor_out_path() -> Path:
    """Support the local showcase and repository robustness-output layouts."""
    candidates = (
        ROOT / "results_showcase" / "submission_robustness" / "leave_one_receptor_out.csv",
        ROOT / "results" / "robustness" / "leave_one_receptor_out.csv",
    )
    return next(path for path in candidates if path.exists())


def save(fig: plt.Figure, name: str) -> None:
    for output in OUTS:
        output.mkdir(parents=True, exist_ok=True)
        fig.savefig(output / f"{name}.png", dpi=350, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def egfr_heatmap() -> None:
    """Render the four prespecified primary-EGFR exclusion effects."""
    labels = ["1M17", "1XKK", "4HJO", "5CAV"]
    values = np.array([[4.24, 3.50, 4.42, 4.79]])
    fig, ax = plt.subplots(figsize=(7.4, 2.8))
    image = ax.imshow(values, cmap="Greens", norm=Normalize(vmin=0, vmax=5), aspect="auto")
    for column, value in enumerate(values[0]):
        ax.text(column, 0, f"{value:+.2f}", ha="center", va="center",
                fontsize=11, fontweight="bold", color="#12321d")
    ax.set_xticks(np.arange(len(labels)), labels, rotation=12)
    ax.set_yticks([0], ["ΔEF1%"])
    ax.set_title("Leave-one-receptor-out robustness of EGFR early enrichment",
                 fontsize=11, fontweight="bold", pad=10)
    ax.tick_params(axis="both", length=0)
    colorbar = fig.colorbar(image, ax=ax, fraction=0.060, pad=0.04)
    colorbar.set_label("Paired EF1% difference")
    save(fig, "figure5_receptor_sensitivity")


def cdk2_forest() -> None:
    """Move the four-primary-receptor CDK2 exclusions to the transfer section."""
    frame = pd.read_csv(leave_one_receptor_out_path())
    receptor_order = ["1fin_a_atp", "2a4l_a_rrc", "1aq1_a_stu", "1pxn_a_ck6"]
    data = frame[(frame.target == "CDK2") & frame.excluded_receptor.isin(receptor_order)].copy()
    data["order"] = pd.Categorical(data.excluded_receptor, receptor_order, ordered=True)
    data = data.sort_values("order")
    labels = [value.split("_")[0].upper() for value in data.excluded_receptor]
    values = data.delta_ef1.to_numpy()
    lower = values - data.delta_ci_lo.to_numpy()
    upper = data.delta_ci_hi.to_numpy() - values
    fig, ax = plt.subplots(figsize=(6.8, 3.0))
    positions = np.arange(len(data))
    ax.errorbar(values, positions, xerr=[lower, upper], fmt="o", color="#296b3d",
                ecolor="#27343a", capsize=3, lw=1.3, ms=6)
    ax.axvline(0, color="#4a565c", lw=1, ls="--")
    ax.set_yticks(positions, labels)
    ax.invert_yaxis()
    ax.set_xlabel("Paired EF1% difference (95% CI)")
    ax.set_title("CDK2 receptor-exclusion sensitivity", fontsize=11, fontweight="bold")
    ax.grid(axis="x", color="#d9e0e3", lw=0.7)
    ax.spines[["top", "right"]].set_visible(False)
    save(fig, "figure6_cdk2_receptor_sensitivity")


if __name__ == "__main__":
    egfr_heatmap()
    cdk2_forest()
