from __future__ import annotations

import html
import json
from pathlib import Path

from egfr_dockingforge.stage12.candidate_card_builder import build_candidate_cards
from egfr_dockingforge.stage12.evidence_aggregation import load_selection_and_inputs
from egfr_dockingforge.stage12.nonclaim_generator import assert_nonclaims, nonclaim_text


def _section(title: str, body: str) -> str:
    return f"<h2>{html.escape(title)}</h2>{body}"


def render_candidate_dossiers(config_path: str | Path) -> dict[str, int]:
    _, paths, _, selection = load_selection_and_inputs(config_path)
    if not any(paths["cards"].glob("*.json")):
        build_candidate_cards(config_path)
    count = 0
    for row in selection[selection["selected_for_detailed_dossier"]].to_dict("records"):
        card = json.loads((paths["cards"] / f"{row['final_candidate_id']}.json").read_text(encoding="utf-8"))
        figure = Path("figures/stage12") / f"{row['final_candidate_id']}_2d.svg"
        body = f"""<!doctype html>
<html><head><meta charset="utf-8"><title>{html.escape(row['final_candidate_id'])}</title>
<style>body{{font-family:Arial,sans-serif;margin:32px;line-height:1.45;color:#172033}}table{{border-collapse:collapse}}td,th{{border:1px solid #cbd5e1;padding:5px 8px}}.note{{background:#f8fafc;border-left:4px solid #64748b;padding:10px}}</style></head>
<body>
<h1>{html.escape(row['final_candidate_id'])}</h1>
<p><strong>Source:</strong> {html.escape(str(row['source']))} | <strong>Novelty:</strong> {html.escape(str(row['novelty_bucket']))} | <strong>Decision:</strong> {html.escape(str(row['decision_label']))}</p>
<p class="note">{html.escape(nonclaim_text())}</p>
{_section("2D Chemistry", f"<img src='../../../{figure.as_posix()}' alt='2D structure'><p><code>{html.escape(str(row['standard_smiles']))}</code></p>")}
{_section("3D Binding Pose", f"<p>Receptor state: {html.escape(str(row['best_receptor_state']))}. Pose file: {html.escape(str(card['best_pose']['pose_file'] or card['best_pose']['missing_pose_reason']))}.</p>")}
{_section("Scores", f"<table><tr><th>Docking</th><th>CNNscore</th><th>CNNaffinity</th><th>Pose confidence</th><th>Calibrated</th></tr><tr><td>{row['best_docking_score']}</td><td>{row['best_gnina_cnnscore']}</td><td>{row['best_gnina_cnnaffinity']}</td><td>{row['best_pose_confidence']}</td><td>{row['best_calibrated_confidence']}</td></tr></table>")}
{_section("Interactions", f"<p>Consensus recall: {row['best_key_interaction_recall_consensus']}; IFP Tanimoto: {row['best_ifp_tanimoto_to_consensus']}.</p>")}
{_section("MD", f"<p>MD status: {html.escape(card['md']['status'])}. Pose stability label: {html.escape(str(card['md']['pose_stability_label']))}. MD not available means no production trajectory was analyzed.</p>")}
{_section("Similarity And Novelty", f"<p>Closest known EGFR ligand: {html.escape(str(row['closest_known_molecule_id']))}; Tanimoto: {row['tanimoto_to_closest_known']}.</p>")}
{_section("Risk", f"<p>{html.escape('; '.join(card['evidence_summary']['main_risks']) or 'No recorded medchem flags.')}</p>")}
{_section("Selection Rationale", f"<p>{html.escape(str(row['selection_reason']))}</p>")}
{_section("Non-Overclaim Statement", f"<p>{html.escape(nonclaim_text())}</p>")}
</body></html>
"""
        assert_nonclaims(body)
        (paths["html"] / f"{row['final_candidate_id']}.html").write_text(body, encoding="utf-8")
        count += 1
    return {"dossiers": count}
