from __future__ import annotations

from pathlib import Path

from egfr_dockingforge.stage12.evidence_aggregation import load_selection_and_inputs
from egfr_dockingforge.stage12.nonclaim_generator import nonclaim_text


def build_dataset_card(config_path: str | Path) -> dict[str, str]:
    _, paths, inputs, _ = load_selection_and_inputs(config_path)
    counts = {name: (0 if frame is None else len(frame)) for name, frame in inputs.items()}
    text = f"""# EGFR DockingForge Project Dataset Card

## Dataset Composition
Stage 1 structural benchmark, Stage 7 candidate library, Stage 8 screening outputs, Stage 9 analog outputs, Stage 10 ablation benchmark outputs, and Stage 11 MD setup/parameterization outputs.

## Source Versions
Input row counts: {counts}.

## Intended Use
Computational EGFR ATP-site docking, pose-confidence triage, interaction-recovery analysis, analog benchmark auditing, and final dossier generation.

## Out-of-Scope Use
Experimental activity, selectivity, toxicity, PK, or clinical inference.

## Provenance
Each Stage 12 card and the reproducibility manifest reference input tables, config hashes, source registry entries, report files, and generation timestamp.

## Licenses And Terms Notes
Upstream structural and ligand sources may carry their own terms. This card records project provenance but does not replace source-specific license review.

## Excluded Records
Records failing project scope, structure quality, candidate filtering, or analog acceptance rules are not promoted as final prospective candidates.

## Source Bias
Current screened candidates are dominated by known-control rows and a small generated analog test set.

## Known Missing Data
Prospective vendor/generated accepted candidates are absent in the current run. Stage 11 production trajectory evidence is not available.

## Limitations
The dataset is designed for computational prioritization and audit reporting only.

{nonclaim_text()}
"""
    path = paths["dataset_cards"] / "project_dataset_card.md"
    path.write_text(text, encoding="utf-8")
    return {"dataset_card": str(path)}
