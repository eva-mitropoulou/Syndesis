from __future__ import annotations

from pathlib import Path

from egfr_dockingforge.stage12.evidence_aggregation import load_selection_and_inputs
from egfr_dockingforge.stage12.nonclaim_generator import nonclaim_text


def build_model_card(config_path: str | Path) -> dict[str, str]:
    _, paths, inputs, _ = load_selection_and_inputs(config_path)
    metrics = inputs.get("stage6_metrics")
    calibration = inputs.get("stage6_calibration")
    text = f"""# Pose Confidence Model Card

## Model Name
EGFR DockingForge Stage 6 pose-confidence model.

## Model Type
Supervised pose-quality confidence model selected by the Stage 6 workflow.

## Training Data
Stage 3-5 redocking/crossdocking pose features and interaction-recovery labels.

## Validation Split
See Stage 6 split artifacts. Metrics rows available in this run: {0 if metrics is None else len(metrics)}.

## Target Label Definition
Pose labels encode computational pose plausibility and interaction recovery, not experimental activity.

## Allowed Use
Ranking and triaging computational docking poses in this EGFR ATP-site workflow.

## Out-of-Scope Use
Experimental activity prediction, cellular potency prediction, selectivity claims, PK, toxicity, or clinical inference.

## Performance Summary
Stage 6 metrics table rows: {0 if metrics is None else len(metrics)}.

## Calibration Summary
Stage 6 calibration table rows: {0 if calibration is None else len(calibration)}.

## Feature Groups
Docking scores, GNINA scores, pose sanity checks, and ProLIF-derived interaction features.

## Leakage Audit Summary
Stage 6 includes leakage-audit checks to avoid directly encoding target labels as model features.

## Limitations
Small EGFR-focused structural benchmark, no wet-lab labels for prospective candidates, and no guarantee of transfer outside the curated receptor ensemble.

## Known Failure Modes
High learned score with weak interaction recovery; false confidence on chemically strained poses; source bias from known-control dominated data.

## Downstream Use
Used by Stage 8, Stage 9 mini-screening, and Stage 12 dossier evidence fields.

{nonclaim_text()}
"""
    path = paths["model_cards"] / "pose_confidence_model_card.md"
    path.write_text(text, encoding="utf-8")
    return {"model_card": str(path)}
