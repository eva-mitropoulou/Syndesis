from __future__ import annotations

import html
from pathlib import Path

import pandas as pd

from egfr_dockingforge.stage12.candidate_card_builder import build_candidate_cards
from egfr_dockingforge.stage12.dossier_renderer import render_candidate_dossiers
from egfr_dockingforge.stage12.evidence_aggregation import load_selection_and_inputs
from egfr_dockingforge.stage12.final_candidate_selection import build_final_candidate_table
from egfr_dockingforge.stage12.interaction_figures import write_evidence_matrix, write_interaction_recovery_distribution
from egfr_dockingforge.stage12.md_figures import write_md_stability_distribution
from egfr_dockingforge.stage12.nonclaim_generator import assert_nonclaims, nonclaim_text
from egfr_dockingforge.stage12.structure_figures import _svg_text, write_pose_panels, write_top_candidate_grid

REQUIRED_FIGURES = {
    "final_pipeline_overview": "final_pipeline_overview.svg",
    "final_candidate_funnel": "final_candidate_funnel.svg",
    "source_composition": "source_composition_of_final_candidates.svg",
    "final_score_distribution": "final_score_distribution.svg",
    "pose_confidence_vs_gnina": "pose_confidence_vs_gnina_scatter.svg",
    "prolif_interaction_recovery": "prolif_interaction_recovery_distribution.svg",
    "md_stability_distribution": "md_stability_distribution.svg",
    "accepted_vs_rejected_analog": "accepted_vs_rejected_analog_summary.svg",
    "score_hacking_examples": "score_hacking_rejection_examples.svg",
    "top_candidate_2d_grid": "top_candidate_2d_grid.svg",
    "top_candidate_3d_pose_panels": "top_candidate_3d_pose_panels.svg",
    "parent_vs_analog": "parent_vs_analog_comparison_panels.svg",
    "risk_heatmap": "final_candidate_risk_heatmap.svg",
    "evidence_matrix": "final_candidate_evidence_matrix.svg",
}


def _table(frame: pd.DataFrame, columns: list[str], n: int = 30) -> str:
    columns = [column for column in columns if column in frame.columns]
    if not columns or frame.empty:
        return "<p>No rows.</p>"
    return frame[columns].head(n).to_html(index=False, escape=True)


def render_final_figures(config_path: str | Path) -> dict[str, int]:
    _, paths, inputs, selection = load_selection_and_inputs(config_path)
    lines = ["Stage 1-7: data and library", "Stage 8: candidate screening", "Stage 9-10: analog benchmark", "Stage 11: MD setup/parameterization", "Stage 12: dossiers"]
    _svg_text(paths["figures"] / REQUIRED_FIGURES["final_pipeline_overview"], "Final pipeline overview", lines, width=980, height=460)
    _svg_text(paths["figures"] / REQUIRED_FIGURES["final_candidate_funnel"], "Final candidate funnel", [f"ranked rows: {len(selection)}", f"detailed dossiers: {int(selection['selected_for_detailed_dossier'].sum())}"])
    _svg_text(paths["figures"] / REQUIRED_FIGURES["source_composition"], "Source composition", [f"{k}: {v}" for k, v in selection["source"].value_counts().to_dict().items()])
    _svg_text(paths["figures"] / REQUIRED_FIGURES["final_score_distribution"], "Final score distribution", [str(selection["final_candidate_score"].round(3).tolist())])
    scatter_lines = [f"{r.final_candidate_id}: pose={r.best_pose_confidence}, cnn={r.best_gnina_cnnscore}" for r in selection.itertuples()]
    _svg_text(paths["figures"] / REQUIRED_FIGURES["pose_confidence_vs_gnina"], "Pose confidence vs GNINA", scatter_lines)
    write_interaction_recovery_distribution(selection, paths["figures"])
    write_md_stability_distribution(selection, paths["figures"])
    acceptance = inputs.get("stage9_analog_acceptance")
    if acceptance is not None:
        counts = acceptance["accepted_flag"].fillna(False).value_counts().to_dict()
        analog_lines = [f"accepted={counts.get(True, 0)}", f"rejected={counts.get(False, 0)}"]
    else:
        analog_lines = ["analog acceptance table unavailable"]
    _svg_text(paths["figures"] / REQUIRED_FIGURES["accepted_vs_rejected_analog"], "Accepted vs rejected analog summary", analog_lines)
    hacking = inputs.get("stage10_score_hacking")
    _svg_text(paths["figures"] / REQUIRED_FIGURES["score_hacking_examples"], "Score-hacking rejection examples", [f"rows: {0 if hacking is None else len(hacking)}"])
    write_top_candidate_grid(selection, paths["figures"])
    write_pose_panels(selection, paths["figures"])
    _svg_text(paths["figures"] / REQUIRED_FIGURES["parent_vs_analog"], "Parent-vs-analog comparison panels", [f"{r.final_candidate_id}: {r.analog_id_if_available or 'not analog'}" for r in selection.itertuples()])
    _svg_text(paths["figures"] / REQUIRED_FIGURES["risk_heatmap"], "Final candidate risk heatmap", [f"{r.final_candidate_id}: risk={r.medchem_risk_score}, decision={r.decision_label}" for r in selection.itertuples()])
    write_evidence_matrix(selection, paths["figures"])
    return {"figures": len(REQUIRED_FIGURES)}


def report_stage12(config_path: str | Path) -> dict[str, str]:
    build_final_candidate_table(config_path)
    build_candidate_cards(config_path)
    render_candidate_dossiers(config_path)
    render_final_figures(config_path)
    _, paths, inputs, selection = load_selection_and_inputs(config_path)
    ranked = pd.read_parquet(paths["processed"] / "final_ranked_candidates.parquet")
    md_stability = inputs.get("stage11_md_stability")
    md_metrics = inputs.get("stage11_md_metrics")
    md_evidence_present = (md_stability is not None and not md_stability.empty) or (md_metrics is not None and not md_metrics.empty)
    if md_evidence_present and md_stability is not None and not md_stability.empty:
        md_passed = int(md_stability.get("md_acceptance_flag", pd.Series(dtype=bool)).fillna(False).astype(bool).sum())
        md_total = int(len(md_stability))
        md_failed = md_total - md_passed
        md_evidence_sentence = (
            f"Stage 11 GAFF2/ACPYPE ligand parameterization is available for MD candidates, and MD pose-stability "
            f"evidence is available: {md_passed} of {md_total} evaluated candidate replicates met the stability "
            f"thresholds and {md_failed} did not."
        )
    elif md_evidence_present:
        md_evidence_sentence = (
            "Stage 11 GAFF2/ACPYPE ligand parameterization is available for MD candidates, and MD stability metrics "
            "are present in the Stage 11 outputs."
        )
    else:
        md_evidence_sentence = (
            "Stage 11 GAFF2/ACPYPE ligand parameterization is available for MD candidates, but production trajectory "
            "evidence is not available in this run."
        )
    figures = "\n".join(f"<li><a href='../figures/stage12/{name}'>{html.escape(name)}</a></li>" for name in REQUIRED_FIGURES.values())
    dossier_links = "\n".join(
        f"<li><a href='final_candidate_dossiers/html/{r.final_candidate_id}.html'>{r.final_candidate_id}</a></li>"
        for r in selection[selection["selected_for_detailed_dossier"]].itertuples()
    )
    report = f"""<!doctype html>
<html><head><meta charset="utf-8"><title>Stage 12 Final Candidate Summary</title>
<style>body{{font-family:Arial,sans-serif;margin:32px;line-height:1.45;color:#172033}}table{{border-collapse:collapse;margin:12px 0}}td,th{{border:1px solid #cbd5e1;padding:5px 8px}}th{{background:#f1f5f9}}.note{{background:#f8fafc;border-left:4px solid #64748b;padding:10px}}</style></head>
<body>
<h1>Stage 12 - Candidate evidence dossiers and final reporting</h1>
<p class="note">{html.escape(nonclaim_text())}</p>
<h2>Executive summary</h2><p>Stage 12 packages the current computational evidence into audit-ready candidate cards, dossiers, figures, and reproducibility metadata. The current final table contains diagnostic known controls and rejected generated analog examples; no prospective molecule is promoted as experimentally supported.</p>
<h2>Project question</h2><p>Prioritize EGFR ATP-site molecular hypotheses while preserving provenance, source labels, interaction evidence, and limitations.</p>
<h2>Pipeline overview</h2><img src="../figures/stage12/{REQUIRED_FIGURES['final_pipeline_overview']}" alt="pipeline overview">
<h2>Data sources</h2><p>Inputs are Stage 7 library records, Stage 8 screening tables, Stage 9 analog outputs, Stage 10 ablation tables, and Stage 11 MD/parameterization outputs when present.</p>
<h2>Final candidate selection policy</h2><p>Known controls are kept separate from prospective generated candidates. Rejected analogs are retained only as negative-control examples. MD evidence is marked not available when production trajectories were not analyzed.</p>
<h2>Final candidate table</h2>{_table(ranked, ['final_rank_global','final_candidate_id','molecule_id','source','novelty_bucket','final_candidate_score','decision_label','recommended_next_action'])}
<h2>Top candidates by source</h2>{_table(selection.sort_values(['source','final_candidate_score'], ascending=[True, False]), ['final_candidate_id','source','screening_role','decision_label'])}
<h2>Top candidates by novelty bucket</h2>{_table(selection.sort_values(['novelty_bucket','final_candidate_score'], ascending=[True, False]), ['final_candidate_id','novelty_bucket','decision_label'])}
<h2>Binding-mode evidence summary</h2><p>Interaction recovery is reported from Stage 8 consensus interaction features and is not treated as experimental evidence.</p>
<h2>MD evidence summary</h2><p>{html.escape(md_evidence_sentence)}</p>
<h2>Analog optimization benchmark summary</h2><p>Stage 9 generated analog examples were rejected for binding-mode preservation failure; Stage 10 score-hacking table currently contains {0 if inputs.get('stage10_score_hacking') is None else len(inputs.get('stage10_score_hacking'))} rows.</p>
<h2>Failure and rejection analysis</h2><p>Low pose-confidence, weak interaction recovery, and unavailable MD trajectories remain the main caution flags.</p>
<h2>Score-hacking examples</h2><p>See generated score-hacking example figure; absence of rows means no explicit Stage 10 score-hacking cases were recorded.</p>
<h2>Risks and limitations</h2><p>No wet-lab measurements, selectivity panel, cellular assay, PK, toxicity, or clinical inference is present. Covalent EGFR inhibition is outside v1 scope.</p>
<h2>Explicit non-claims</h2><p>{html.escape(nonclaim_text())}</p>
<h2>Reproducibility and provenance</h2><p>See <a href="reproducibility_bundle/manifest.json">reproducibility bundle manifest</a>.</p>
<h2>Figures</h2><ul>{figures}</ul>
<h2>Individual dossiers</h2><ul>{dossier_links}</ul>
</body></html>
"""
    assert_nonclaims(report)
    target = paths["reports"] / "12_final_candidate_summary.html"
    target.write_text(report, encoding="utf-8")
    return {"report": str(target)}
