#!/usr/bin/env python3
"""Verify the principal paper values from frozen Syndesis result tables."""

from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]


def require_close(value: float, expected: float, label: str, tolerance: float = 1e-6) -> None:
    if abs(value - expected) > tolerance:
        raise SystemExit(f"{label}: expected {expected}, found {value}")


def main() -> None:
    paper = pd.read_csv(ROOT / "results/statistics/paper_metrics.csv")
    egfr_gnina = paper.loc[
        paper["target"].eq("EGFR") & paper["arm"].eq("gnina") & paper["metric"].eq("ef1"),
        "estimate",
    ].iloc[0]
    egfr_coupled = paper.loc[
        paper["target"].eq("EGFR") & paper["arm"].eq("coupled") & paper["metric"].eq("ef1"),
        "estimate",
    ].iloc[0]
    require_close(egfr_gnina, 11.98, "EGFR GNINA EF1%", 0.01)
    require_close(egfr_coupled, 16.40, "EGFR pose-coupled EF1%", 0.01)

    cdk2 = pd.read_csv(ROOT / "results/statistics/cdk2/top1_active_counts.csv")
    cdk2_gnina = cdk2.loc[cdk2["arm"].eq("gnina"), "ef1"].iloc[0]
    cdk2_coupled = cdk2.loc[cdk2["arm"].eq("coupled"), "ef1"].iloc[0]
    require_close(cdk2_gnina, 11.6017, "CDK2 GNINA EF1%", 1e-4)
    require_close(cdk2_coupled, 14.9768, "CDK2 pose-coupled EF1%", 1e-4)

    effects = pd.read_csv(ROOT / "results/statistics/cdk2/paired_metric_effects.csv")
    effect = effects.loc[
        effects["metric"].eq("ef1") & effects["contrast"].eq("coupled_minus_gnina"),
        "estimate",
    ].iloc[0]
    require_close(effect, 3.3750503198103505, "CDK2 paired EF1% effect")
    print("Verified frozen EGFR and corrected CDK2 paper values.")


if __name__ == "__main__":
    main()
