from __future__ import annotations

from pathlib import Path

import pandas as pd


def _bar_svg(labels: list[str], values: list[float], title: str) -> str:
    width, height = 900, 320
    maxv = max(values + [1e-9])
    bars = []
    step = width / max(len(values), 1)
    for i, (label, value) in enumerate(zip(labels, values)):
        h = 220 * value / maxv if maxv else 0
        x = i * step + 10
        y = 260 - h
        bars.append(f'<rect x="{x:.1f}" y="{y:.1f}" width="{max(step-20,8):.1f}" height="{h:.1f}" fill="#4477aa"/><text x="{x:.1f}" y="285" font-size="10" transform="rotate(30 {x:.1f},285)">{label}</text><text x="{x:.1f}" y="{y-4:.1f}" font-size="10">{value:.2f}</text>')
    return f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}"><text x="20" y="25" font-size="18">{title}</text>{"".join(bars)}</svg>'


def make_ablation_plots(strategy_metrics: pd.DataFrame, paths: dict[str, Path]) -> list[Path]:
    specs = [
        ("accepted_rate.svg", "accepted_analog_rate_pre_md", "Accepted analog rate"),
        ("score_hacking_rate.svg", "score_hacking_rate", "Score-hacking rate"),
        ("validity_rate.svg", "validity_rate", "Validity rate"),
        ("delta_pose_confidence.svg", "mean_delta_pose_confidence", "Mean delta pose confidence"),
        ("delta_key_interaction.svg", "mean_delta_key_interaction_recall", "Mean delta key interaction recall"),
        ("pareto_accepted_vs_hacking.svg", "accepted_analog_rate_pre_md", "Accepted rate Pareto proxy"),
        ("pareto_accepted_vs_cost.svg", "accepted_analog_rate_pre_md", "Accepted rate vs compute proxy"),
        ("ablation_waterfall.svg", "accepted_analog_rate_pre_md", "Ablation waterfall proxy"),
        ("seed_paired_comparison.svg", "accepted_analog_rate_pre_md", "Seed paired comparison proxy"),
        ("rejection_reasons.svg", "binding_mode_break_rate", "Binding-mode break/rejection rate"),
    ]
    out = []
    labels = strategy_metrics["strategy_name"].tolist()
    for filename, col, title in specs:
        path = paths["figures"] / filename
        vals = strategy_metrics[col].fillna(0).astype(float).tolist()
        path.write_text(_bar_svg(labels, vals, title), encoding="utf-8")
        out.append(path)
    return out
