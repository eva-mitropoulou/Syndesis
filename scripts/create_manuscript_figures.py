#!/usr/bin/env python3
"""Create publication figures for the Syndesis manuscript from tracked showcase data."""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import TwoSlopeNorm
from rdkit import Chem
from rdkit.Chem import Draw


ROOT = Path(__file__).resolve().parents[1]
SHOWCASE = ROOT / "results_showcase"
OUT = ROOT / "figures" / "manuscript"
OUT.mkdir(parents=True, exist_ok=True)

COLORS = {
    "dock": "#637381",
    "inter": "#2596a8",
    "gnina": "#e07a32",
    "combined": "#296b3d",
    "reject": "#a33a3a",
    "stable": "#296b3d",
}
plt.rcParams.update(
    {
        "font.family": "DejaVu Sans",
        "font.size": 9,
        "axes.titleweight": "bold",
        "axes.spines.top": False,
        "axes.spines.right": False,
        "figure.dpi": 160,
        "savefig.dpi": 350,
    }
)


def save(fig, name: str) -> None:
    fig.savefig(OUT / f"{name}.png", bbox_inches="tight", facecolor="white")
    fig.savefig(OUT / f"{name}.pdf", bbox_inches="tight", facecolor="white")
    plt.close(fig)


def enrichment() -> None:
    metrics_table = pd.read_csv(SHOWCASE / "submission_robustness" / "bootstrap_metric_intervals.csv")
    targets = ["EGFR", "CDK2"]
    arms = ["gnina", "coupled"]
    labels = ["GNINA", "GNINA +\ninteraction"]
    colors = [COLORS["gnina"], COLORS["combined"]]
    metrics = [("roc_auc", "ROC-AUC", (0.5, 0.85)), ("ef1", "EF1%", (0, 19)), ("bedroc_80_5", "BEDROC", (0, 0.34))]
    fig, axes = plt.subplots(1, 3, figsize=(10.6, 3.4), sharex=True)
    positions = np.array([0, 1])
    width = 0.18
    for ax, (column, title, ylim) in zip(axes, metrics):
        for arm_index, (arm, label, color) in enumerate(zip(arms, labels, colors)):
            values, lower, upper = [], [], []
            for target in targets:
                row = metrics_table[(metrics_table.target == target) & (metrics_table.arm == arm) & (metrics_table.metric == column)].iloc[0]
                values.append(row.estimate)
                lower.append(row.estimate - row.ci_lo)
                upper.append(row.ci_hi - row.estimate)
            offset = (arm_index - 0.5) * width
            ax.bar(positions + offset, values, width, color=color, label=label, zorder=2)
            ax.errorbar(positions + offset, values, yerr=[lower, upper], fmt="none", color="#27343a", capsize=2, lw=0.8, zorder=3)
        ax.set_title(title)
        ax.set_xticks(positions, targets)
        ax.set_ylim(*ylim)
        ax.grid(axis="y", color="#d9e0e3", lw=0.7, zorder=0)
    handles, legend_labels = axes[-1].get_legend_handles_labels()
    fig.legend(handles, legend_labels, loc="lower center", ncol=2, frameon=False, bbox_to_anchor=(0.5, -0.12))
    fig.suptitle("Pose-coupled interaction weighting improves EGFR early enrichment", y=1.03, fontsize=11, fontweight="bold")
    save(fig, "figure2_enrichment")


def paired_deltas() -> None:
    final_dir = SHOWCASE / "submission_robustness" / "cdk2_four_receptor_four_prior"
    deltas = pd.read_csv(final_dir / "paired_metric_effects.csv")
    formula = pd.read_csv(final_dir / "interaction_formula_sensitivity.csv")
    cdk2_prior = pd.read_csv(final_dir / "cdk2_native_prior_sensitivity.csv")
    cdk2_gnina = pd.read_csv(
        SHOWCASE / "submission_robustness" / "cdk2_late_fusion_metrics.csv"
    )

    # Display observed full-dataset differences. Paired-bootstrap draws provide
    # the uncertainty intervals but their means are not substituted as points.
    egfr_baseline = formula[
        (formula.target == "EGFR")
        & (formula.combination == "multiplicative")
        & (formula.term == "recall")
        & (formula["lambda"] == 0.0)
    ].iloc[0]
    egfr_coupled = formula[
        (formula.target == "EGFR")
        & (formula.combination == "multiplicative")
        & (formula.term == "recall")
        & (formula["lambda"] == 1.0)
    ].iloc[0]
    cdk2_coupled = cdk2_prior[
        cdk2_prior.prior_definition == "primary_union_recall"
    ].iloc[0]
    cdk2_baseline = cdk2_gnina[cdk2_gnina.arm == "gnina"].set_index("metric")

    observed = {
        ("EGFR", "ef1"): egfr_coupled.ef1 - egfr_baseline.ef1,
        ("EGFR", "bedroc_80_5"): egfr_coupled.bedroc - egfr_baseline.bedroc,
        ("CDK2", "ef1"): cdk2_coupled.ef1 - cdk2_baseline.loc["ef1", "estimate"],
        ("CDK2", "bedroc_80_5"): (
            cdk2_coupled.bedroc - cdk2_baseline.loc["bedroc_80_5", "estimate"]
        ),
    }
    fig, axes = plt.subplots(1, 2, figsize=(8.4, 3.6))
    panels = [("ef1", "EF1% difference"), ("bedroc", "BEDROC difference")]
    for ax, (metric, title) in zip(axes, panels):
        rows = []
        for target in ("EGFR", "CDK2"):
            key = "bedroc_80_5" if metric == "bedroc" else metric
            row = deltas[(deltas.target == target) & (deltas.metric == key)].iloc[0]
            # Use the observed full-dataset effect with the paired-bootstrap CI.
            point = observed[(target, key)]
            low = row.ci_lo
            high = row.ci_hi
            if target == "CDK2" and key == "ef1":
                low = cdk2_coupled.delta_ef1_ci_lo
                high = cdk2_coupled.delta_ef1_ci_hi
            rows.append((target, point, low, high))
        points = [row[1] for row in rows]
        lows = [row[2] for row in rows]
        highs = [row[3] for row in rows]
        span = max(highs) - min(lows)
        padding = max(span * 0.12, 0.012 if metric == "bedroc" else 0.20)
        ax.set_xlim(min(lows) - padding, max(highs) + padding)
        for y, (target, point, low, high) in enumerate(rows[::-1]):
            ax.errorbar(point, y, xerr=[[point - low], [high - point]], fmt="o", color=COLORS["combined"], capsize=4, lw=1.5)
            ax.annotate(
                f"{point:+.3f}" if metric == "bedroc" else f"{point:+.2f}",
                xy=(point, y),
                xytext=(7, 7),
                textcoords="offset points",
                ha="left",
                va="bottom",
                fontsize=8,
                color="#172126",
            )
        ax.axvline(0, color="#4a565c", lw=1, ls="--")
        ax.set_yticks(range(len(rows)), [row[0] for row in rows[::-1]])
        ax.set_title(title)
        ax.grid(axis="x", color="#d9e0e3", lw=0.7)
        ax.set_xlabel("GNINA + interaction minus GNINA (95% CI)")
    fig.suptitle("Paired-bootstrap effects of pose-coupled interaction weighting", y=1.03, fontsize=11, fontweight="bold")
    save(fig, "figure3_paired_deltas")


def permutation_nulls() -> None:
    summary = pd.read_csv(SHOWCASE / "submission_robustness" / "permutation_null_summary.csv")
    draws = pd.read_parquet(SHOWCASE / "submission_robustness" / "permutation_null_draws.parquet")
    egfr_summary = pd.read_csv(SHOWCASE / "submission_robustness" / "egfr_atp_prior_permutation_summary.csv")
    egfr_draws = pd.read_parquet(SHOWCASE / "submission_robustness" / "egfr_atp_prior_permutation_draws.parquet")
    labels = {
        "all_ligand": "All-ligand",
        "heavy_atom_decile": "Heavy-atom-count-matched",
        "class_conditional": "Class-conditional",
    }
    colors = {"all_ligand": "#637381", "heavy_atom_decile": "#2596a8", "class_conditional": "#e07a32"}
    fig, axes = plt.subplots(1, 2, figsize=(10.2, 3.8), sharey=False)
    for ax, target in zip(axes, ["EGFR", "CDK2"]):
        target_summary = egfr_summary if target == "EGFR" else summary[summary.target == target]
        target_draws = egfr_draws if target == "EGFR" else draws[draws.target == target]
        for null_name in labels:
            values = target_draws[target_draws["null"] == null_name]["ef1"]
            ax.hist(values, bins=28, density=True, histtype="step", lw=1.6,
                    color=colors[null_name], label=labels[null_name])
        observed = target_summary.iloc[0].observed_ef1
        ax.axvline(observed, color=COLORS["combined"], lw=2.0, label="Observed" if target == "EGFR" else None)
        ax.set_title(target)
        ax.set_xlabel("Permuted EF1%")
        ax.set_ylabel("Density")
        ax.grid(axis="y", color="#d9e0e3", lw=0.7)
    handles, legend_labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, legend_labels, ncol=4, frameon=False, loc="lower center", bbox_to_anchor=(0.5, -0.08))
    fig.suptitle("Observed early enrichment exceeds three interaction-assignment nulls", y=1.02, fontsize=11, fontweight="bold")
    save(fig, "figure4_permutation_nulls")


def receptor_sensitivity() -> None:
    frame = pd.read_csv(SHOWCASE / "submission_robustness" / "leave_one_receptor_out.csv")
    egfr = pd.read_csv(SHOWCASE / "submission_robustness" / "egfr_atp_prior_leave_one_receptor_out.csv")
    targets = ["EGFR", "CDK2"]
    fig, axes = plt.subplots(1, 2, figsize=(10.2, 4.1))
    for ax, target in zip(axes, targets):
        data = egfr.sort_values("excluded_receptor") if target == "EGFR" else frame[frame.target == target].sort_values("excluded_receptor")
        labels = [value.split("_")[0].upper() for value in data.excluded_receptor]
        values = data.delta_ef1.to_numpy()[None, :]
        image = ax.imshow(values, cmap="RdYlGn", norm=TwoSlopeNorm(vmin=-max(5, abs(values.min())), vcenter=0, vmax=max(5, abs(values.max()))), aspect="auto")
        for column, row in enumerate(data.itertuples()):
            ax.text(column, 0, f"{row.delta_ef1:+.2f}", ha="center", va="center",
                    fontsize=9, fontweight="bold", color="#172126")
        ax.set_xticks(np.arange(len(labels)), labels)
        ax.set_yticks([0], ["ΔEF1%"] if target == "EGFR" else [""])
        ax.set_title(target)
        ax.tick_params(axis="x", rotation=35)
        fig.colorbar(image, ax=ax, fraction=0.035, pad=0.03, label="EF1% difference")
    fig.suptitle("Leave-one-receptor-out sensitivity", y=1.02, fontsize=11, fontweight="bold")
    save(fig, "figure5_receptor_sensitivity")


def md_stability() -> None:
    summary = pd.read_parquet(ROOT / "data/processed/stage11/md_candidate_summary.parquet")
    frame = summary[summary.num_replicates_completed > 0].copy()
    names = {
        "mdcand_001": "Control 001", "mdcand_002": "Control 002", "mdcand_003": "Control 003",
        "mdcand_004": "Analog 004", "mdcand_005": "Analog 005", "mdcand_006": "Analog 006",
        "mdcand_neg01": "Mis-docked control",
    }
    frame["finalist"] = frame.md_candidate_id.map(names)
    frame["rmsd"] = frame.median_ligand_rmsd
    frame["key"] = frame.median_key_interaction_occupancy
    frame["verdict"] = np.where(frame.final_md_label.eq("md_stable"), "stable", "unstable")
    metrics = pd.read_parquet(ROOT / "data/processed/stage11/md_metrics.parquet")
    spread = metrics.groupby("md_candidate_id").ligand_rmsd_median_angstrom.agg(["min", "max"])
    frame = frame.join(spread, on="md_candidate_id")
    fig, axes = plt.subplots(1, 2, figsize=(9.2, 4.1), gridspec_kw={"width_ratios": [1.1, 1]})
    ordered = frame.sort_values("rmsd", ascending=False).reset_index(drop=True)
    ys = np.arange(len(ordered))
    colors = [COLORS["stable"] if item == "stable" else COLORS["reject"] for item in ordered.verdict]
    axes[0].errorbar(ordered.rmsd, ys, xerr=[ordered.rmsd - ordered["min"], ordered["max"] - ordered.rmsd],
                     fmt="none", color="#27343a", capsize=2, zorder=2)
    axes[0].scatter(ordered.rmsd, ys, s=50, c=colors, zorder=3)
    axes[0].axvline(3.0, color="#4a565c", ls="--", lw=1)
    axes[0].set_yticks(ys, ordered.finalist)
    axes[0].set_xlabel("Ligand RMSD in binding-site frame (Å)")
    axes[0].set_title("Pose geometry")
    axes[0].grid(axis="x", color="#d9e0e3", lw=0.7)
    axes[0].text(3.05, 0.03, "3 Å gate", fontsize=7.5, va="bottom", color="#4a565c",
                 transform=axes[0].get_xaxis_transform())
    axes[1].barh(ys, ordered.key, height=0.52, label="Mean key contacts", color="#e07a32")
    axes[1].axvline(0.50, color="#e07a32", ls="--", lw=1)
    axes[1].set_xlim(0, 1.05)
    axes[1].set_yticks(ys, [])
    axes[1].set_xlabel("Interaction occupancy")
    axes[1].set_title("Interaction persistence")
    axes[1].grid(axis="x", color="#d9e0e3", lw=0.7)
    axes[1].text(0.51, 0.03, "0.50 gate", fontsize=7.5, va="bottom", color="#4a565c",
                 transform=axes[1].get_xaxis_transform())
    fig.suptitle("The MD gate separates majority-stable systems from the deliberately mis-docked control in this selected set", y=1.02, fontsize=11, fontweight="bold")
    save(fig, "figure6_md_stability")


def prospective() -> None:
    frame = pd.read_parquet("/mnt/e/egfr_prospective/prospective_ranked_corrected.parquet")
    audit = pd.read_csv(SHOWCASE / "submission_robustness" / "prospective_gate_audit.csv").iloc[0]
    frame["passes_gate"] = ((frame.cnnscore >= audit.cnn_top_decile_threshold) &
                             (frame.key_interaction_recall_consensus >= audit.global_median_recall_threshold))
    accepted = frame[frame.passes_gate].sort_values("gnina_interaction_score", ascending=False)
    top = accepted.head(6).copy()
    fig, axes = plt.subplots(1, 2, figsize=(10.6, 4.5), gridspec_kw={"width_ratios": [1.25, 1]})
    ax = axes[0]
    ax.scatter(frame.cnnscore, frame.key_interaction_recall_consensus, s=10, color="#b9c3c8", alpha=0.45, linewidths=0)
    scatter = ax.scatter(
        accepted.cnnscore,
        accepted.key_interaction_recall_consensus,
        s=24 + 45 * accepted.qed,
        c=accepted.qed,
        cmap="viridis",
        edgecolors="#27343a",
        linewidths=0.35,
        alpha=0.92,
    )
    ax.axvline(audit.cnn_top_decile_threshold, color="#4a565c", ls="--", lw=0.9)
    ax.axhline(audit.global_median_recall_threshold, color="#4a565c", ls="--", lw=0.9)
    ax.set_xlabel("GNINA CNNscore")
    ax.set_ylabel("Native-union interaction recall")
    ax.set_title("Conjunctive gate and coupled ranking")
    colorbar = fig.colorbar(scatter, ax=ax, pad=0.02)
    colorbar.set_label("QED")
    ax.grid(color="#d9e0e3", lw=0.7)
    structure_ax = axes[1]
    structure_ax.set_axis_off()
    molecules = [Chem.MolFromSmiles(value) for value in top.smiles]
    legends = [f"P{i}: {zinc}" for i, zinc in enumerate(top.zinc_id, start=1)]
    grid = Draw.MolsToGridImage(molecules, molsPerRow=2, subImgSize=(300, 220), legends=legends, useSVG=False)
    structure_ax.imshow(grid)
    structure_ax.set_title("Highest-ranked gate-passing structures")
    fig.suptitle("Prospective EGFR prioritization yields testable, unlabeled hypotheses", y=1.02, fontsize=11, fontweight="bold")
    save(fig, "figure7_prospective")


if __name__ == "__main__":
    enrichment()
    paired_deltas()
    permutation_nulls()
    receptor_sensitivity()
    md_stability()
    prospective()
