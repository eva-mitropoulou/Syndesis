from __future__ import annotations

from pathlib import Path

from egfr_dockingforge.stage12.candidate_card_builder import build_candidate_cards
from egfr_dockingforge.stage12.dataset_card import build_dataset_card
from egfr_dockingforge.stage12.dossier_renderer import render_candidate_dossiers
from egfr_dockingforge.stage12.final_candidate_selection import build_final_candidate_table
from egfr_dockingforge.stage12.model_card import build_model_card
from egfr_dockingforge.stage12.provenance_bundle import build_provenance_bundle
from egfr_dockingforge.stage12.report_stage12 import render_final_figures, report_stage12


def run_stage12_all(config_path: str | Path) -> dict[str, str]:
    build_final_candidate_table(config_path)
    build_candidate_cards(config_path)
    render_candidate_dossiers(config_path)
    render_final_figures(config_path)
    build_model_card(config_path)
    build_dataset_card(config_path)
    report = report_stage12(config_path)
    bundle = build_provenance_bundle(config_path)
    return {"report": report["report"], "manifest": bundle["manifest"]}
