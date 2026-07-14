from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

from syndesis.common.io import write_table


def _stage11_md_ablation_row(inputs: dict[str, pd.DataFrame | None] | None) -> dict:
    base = {
        "ablation_id": "md_filter_ablation",
        "base_strategy": "rdkit_rule_based",
        "removed_component": "",
        "added_component": "md_filter",
        "accepted_rate_change": 0.0,
        "score_hacking_rate_change": 0.0,
        "binding_mode_break_rate_change": 0.0,
        "medchem_risk_change": 0.0,
        "runtime_change": 0.0,
    }
    if not inputs:
        return {
            **base,
            "conclusion": "pending_stage11",
            "evidence_strength": "pending",
            "warnings_json": json.dumps(["stage11_md_metrics_missing"]),
        }

    metrics = inputs.get("stage11_md_metrics")
    labels = inputs.get("stage11_md_pose_stability_labels")
    if metrics is None or metrics.empty or labels is None or labels.empty:
        return {
            **base,
            "conclusion": "pending_stage11",
            "evidence_strength": "pending",
            "warnings_json": json.dumps(["stage11_md_metrics_missing"]),
        }

    complete = int(metrics.get("trajectory_analysis_status", pd.Series(dtype=str)).eq("complete").sum())
    accepted = int(labels.get("md_acceptance_flag", pd.Series(dtype=bool)).fillna(False).astype(bool).sum())
    total = int(len(labels))
    if complete < total:
        conclusion = "stage11_md_incomplete"
        evidence_strength = "partial_md_evidence"
        warnings = ["stage11_md_incomplete"]
    elif accepted == 0:
        conclusion = "stage11_completed_no_candidates_passed_md_filter"
        evidence_strength = "md_filter_evidence_available"
        warnings = []
    else:
        conclusion = "stage11_completed_some_candidates_passed_md_filter"
        evidence_strength = "md_filter_evidence_available"
        warnings = []

    # Pre-MD accepted rate must reflect the actual stage 9/10 acceptance BEFORE the
    # stage 11 MD filter, not a fabricated 1.0. Derive it from the stage9 acceptance
    # table (num accepted / num evaluated) when available; if it cannot be derived or
    # is undefined, fall back to 0.0 (never 1.0) and record a warning.
    acceptance = inputs.get("stage9_analog_acceptance")
    if acceptance is not None and not acceptance.empty:
        evaluated = int(len(acceptance))
        accepted_pre_md = int(acceptance.get("accepted_flag", pd.Series(dtype=bool)).fillna(False).astype(bool).sum())
        pre_rate = accepted_pre_md / evaluated if evaluated else 0.0
    else:
        pre_rate = 0.0
        warnings = [*warnings, "pre_md_accept_rate_unavailable_defaulted_to_zero"]
    post_rate = accepted / max(total, 1)
    return {
        **base,
        "accepted_rate_change": post_rate - pre_rate,
        "conclusion": conclusion,
        "evidence_strength": evidence_strength,
        "warnings_json": json.dumps(warnings),
    }


def build_ablation_summary(
    strategy_metrics: pd.DataFrame,
    paths: dict[str, Path],
    inputs: dict[str, pd.DataFrame | None] | None = None,
) -> pd.DataFrame:
    idx = strategy_metrics.set_index("strategy_name")
    pairs = [
        ("gnina_vs_docking", "docking_score_only_optimization", "gnina_only_optimization", "gnina"),
        ("constrained_rdkit_vs_gnina", "gnina_only_optimization", "rdkit_rule_based", "interaction_and_pose_constraints"),
    ]
    rows = []
    for aid, base, comp, added in pairs:
        if base in idx.index and comp in idx.index:
            b = idx.loc[base]
            c = idx.loc[comp]
            rows.append({
                "ablation_id": aid, "base_strategy": base, "removed_component": "", "added_component": added,
                "accepted_rate_change": c["accepted_analog_rate_pre_md"] - b["accepted_analog_rate_pre_md"],
                "score_hacking_rate_change": c["score_hacking_rate"] - b["score_hacking_rate"],
                "binding_mode_break_rate_change": c["binding_mode_break_rate"] - b["binding_mode_break_rate"],
                "medchem_risk_change": c["medchem_risk_rejection_rate"] - b["medchem_risk_rejection_rate"],
                "runtime_change": 0.0,
                "conclusion": "no_difference_observed" if c["accepted_analog_rate_pre_md"] == b["accepted_analog_rate_pre_md"] else "rate_changed",
                "evidence_strength": "exploratory_small_seed_count",
                "warnings_json": json.dumps([]),
            })
    rows.append(_stage11_md_ablation_row(inputs))
    out = pd.DataFrame(rows)
    write_table(paths["processed"] / "ablation_summary.parquet", out)
    write_table(paths["processed"] / "ablation_summary.csv", out)
    return out
