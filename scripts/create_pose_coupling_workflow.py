#!/usr/bin/env python3
"""Render the clean Figure 1 workflow for pose-coupled ranking."""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch


ROOT = Path(__file__).resolve().parents[1]
MANUSCRIPT_OUT = ROOT / "manuscript" / "figures"
ARCHIVE_OUT = ROOT / "figures" / "manuscript"


def box(ax, y, text, color, *, height=0.085, fontsize=15):
    """Draw a workflow stage and return its vertical bounds."""
    ax.add_patch(FancyBboxPatch(
        (0.13, y), 0.74, height,
        boxstyle="round,pad=0.012,rounding_size=0.025",
        facecolor=color, edgecolor="#18324a", linewidth=1.4,
    ))
    ax.text(0.50, y + height / 2, text, ha="center", va="center",
            fontsize=fontsize, color="#10263b", weight="semibold",
            linespacing=1.35)
    return y, y + height


def arrow(ax, upper, lower):
    ax.annotate("", xy=(0.50, lower), xytext=(0.50, upper),
                arrowprops={"arrowstyle": "-|>", "color": "#40627d", "lw": 2.0,
                           "shrinkA": 4, "shrinkB": 5})


def render(output: Path) -> None:
    fig, ax = plt.subplots(figsize=(8.5, 11.0), dpi=300)
    fig.patch.set_facecolor("white")
    ax.set(xlim=(0, 1), ylim=(0, 1))
    ax.axis("off")
    ax.text(0.50, 0.966, "Pose-coupled interaction-aware ranking",
            ha="center", va="center", fontsize=17, weight="bold", color="#10263b")
    ax.text(0.50, 0.932,
            "Interaction evidence contributes only when it belongs to the same receptor-specific pose as CNNscore.",
            ha="center", va="center", fontsize=8.5, color="#496276")

    stages = [
        box(ax, 0.826, "Four EGFR ATP-site native complexes\n1M17/AQ4  •  1XKK/FMM  •  4HJO/AQ4  •  5CAV/4ZQ", "#dbeafe", height=0.087, fontsize=10.5),
        box(ax, 0.698, "Residue–interaction fingerprints", "#e0f2fe", fontsize=12.5),
        box(ax, 0.570, "Union interaction prior  C", "#ccfbf1", fontsize=12.5),
        box(ax, 0.424, "Candidate docked independently into\nfour receptor states", "#e0f2fe", height=0.090, fontsize=12.5),
        box(ax, 0.278, "For each selected pose:\nCNNscore(i,r) and recall R(i,r)", "#fef3c7", height=0.090, fontsize=11.7),
        box(ax, 0.143, "Same-pose coupled score\nCNNscore(i,r) × [1 + R(i,r)]", "#fde68a", height=0.090, fontsize=12.5),
        box(ax, 0.030, "Ligand ranking:\nmaximum coupled score across receptor states", "#bbf7d0", height=0.072, fontsize=10.8),
    ]
    for earlier, later in zip(stages, stages[1:]):
        arrow(ax, earlier[0], later[1])
    ax.text(0.50, 0.118, "CNNscore and recall are never pooled across different poses.",
            ha="center", va="center", fontsize=8.3, color="#7c4a03", style="italic")

    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)


if __name__ == "__main__":
    render(MANUSCRIPT_OUT / "figure1_pose_coupling_workflow.png")
    render(ARCHIVE_OUT / "figure1_pose_coupling_workflow.png")
