from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
import pytest

from syndesis.stage12.candidate_card_builder import build_candidate_cards
from syndesis.stage12.dossier_renderer import render_candidate_dossiers
from syndesis.stage12.final_candidate_selection import build_final_candidate_table
from syndesis.stage12.nonclaim_generator import REQUIRED_NON_CLAIMS


CONFIG = Path("configs/stage12_candidate_dossiers.yaml")
pytestmark = pytest.mark.integration


def test_selected_candidates_have_cards_dossiers_figures_and_pose_status() -> None:
    build_final_candidate_table(CONFIG)
    build_candidate_cards(CONFIG)
    render_candidate_dossiers(CONFIG)
    selection = pd.read_parquet("data/processed/stage12/final_candidate_selection.parquet")
    selected = selection[selection["selected_for_detailed_dossier"]]
    assert not selected.empty
    for row in selected.to_dict("records"):
        card_path = Path("reports/final_candidate_dossiers/cards") / f"{row['final_candidate_id']}.json"
        html_path = Path("reports/final_candidate_dossiers/html") / f"{row['final_candidate_id']}.html"
        figure_path = Path("figures/stage12") / f"{row['final_candidate_id']}_2d.svg"
        assert card_path.exists()
        assert html_path.exists()
        assert figure_path.exists()
        card = json.loads(card_path.read_text(encoding="utf-8"))
        assert card["best_pose"]["pose_file"] or card["best_pose"]["missing_pose_reason"]
        text = html_path.read_text(encoding="utf-8")
        for claim in REQUIRED_NON_CLAIMS:
            assert claim in text
