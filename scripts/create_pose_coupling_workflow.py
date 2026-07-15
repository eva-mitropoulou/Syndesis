#!/usr/bin/env python3
"""Render the clean Figure 1 workflow for pose-coupled ranking."""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch


ROOT = Path(__file__).resolve().parents[1]
MANUSCRIPT_OUT = ROOT / "manuscript" / "figures"
ARCHIVE_OUT = ROOT / "figures" / "manuscript"


def box(ax, y, text, color, *, height=0.075, fontsize=8.2, text_color="#18324a"):
    """Draw a workflow stage and return its vertical bounds."""
    ax.add_patch(FancyBboxPatch(
        (0.13, y), 0.74, height,
        boxstyle="round,pad=0.012,rounding_size=0.025",
        facecolor=color, edgecolor="#40627d", linewidth=1.0, zorder=2,
    ))
    ax.text(0.50, y + height / 2, text, ha="center", va="center",
            fontsize=fontsize, color=text_color, weight="semibold",
            linespacing=1.2, zorder=3)
    return y, y + height


def arrow(ax, upper, lower):
    ax.annotate("", xy=(0.50, lower), xytext=(0.50, upper),
                arrowprops={"arrowstyle": "-|>", "color": "#61788b", "lw": 1.25,
                           "shrinkA": 3, "shrinkB": 3.5}, zorder=1)


def render(output: Path) -> None:
    fig, ax = plt.subplots(figsize=(7.2, 8.4), dpi=300)
    fig.patch.set_facecolor("white")
    ax.set(xlim=(0, 1), ylim=(0, 1))
    ax.axis("off")

    ax.text(0.50, 0.965, "Pose-coupled interaction-aware ranking",
            ha="center", va="center", fontsize=10.5, weight="bold", color="#18324a")

    stages = [
        # Each stage has a distinct muted fill so that colour never implies an
        # unstated equivalence between separate workflow steps.
        box(ax, 0.827, "Four EGFR ATP-site native complexes\n1M17/AQ4  •  1XKK/FMM  •  4HJO/AQ4  •  5CAV/4ZQ", "#e8f0f7", height=0.078, fontsize=7.3),
        box(ax, 0.704, "Residue–interaction fingerprints", "#eeeaf6", fontsize=8.2),
        box(ax, 0.581, "Union interaction prior  C", "#e1f0ea", fontsize=8.2),
        box(ax, 0.438, "Candidate docked independently into\nfour receptor states", "#e5f1f3", height=0.080, fontsize=8.0),
        box(ax, 0.315, "For each selected pose:\nCNNscore(i,r) and recall R(i,r)", "#f3f5f6", height=0.080, fontsize=7.6),
        box(ax, 0.192, "Same-pose coupled score\nCNNscore(i,r) × [1 + R(i,r)]", "#f2e6c9", height=0.080, fontsize=8.0),
        box(ax, 0.055, "Ligand ranking:\nmaximum coupled score across receptor states", "#38666a", height=0.065, fontsize=7.2, text_color="white"),
    ]
    for earlier, later in zip(stages, stages[1:]):
        arrow(ax, earlier[0], later[1])

    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)


if __name__ == "__main__":
    render(MANUSCRIPT_OUT / "figure1_pose_coupling_workflow.png")
    render(ARCHIVE_OUT / "figure1_pose_coupling_workflow.png")
