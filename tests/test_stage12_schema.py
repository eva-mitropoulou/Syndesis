from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
import pytest

from egfr_dockingforge.stage12.candidate_card_schema import validate_candidate_card
from egfr_dockingforge.stage12.candidate_dossiers import run_stage12_all
from egfr_dockingforge.stage12.schemas import FINAL_CANDIDATE_COLUMNS, FINAL_RANKED_COLUMNS


CONFIG = Path("configs/stage12_candidate_dossiers.yaml")
pytestmark = pytest.mark.integration


def test_stage12_outputs_match_required_schemas() -> None:
    run_stage12_all(CONFIG)
    selection = pd.read_parquet("data/processed/stage12/final_candidate_selection.parquet")
    ranked = pd.read_parquet("data/processed/stage12/final_ranked_candidates.parquet")
    assert set(FINAL_CANDIDATE_COLUMNS).issubset(selection.columns)
    assert set(FINAL_RANKED_COLUMNS).issubset(ranked.columns)
    assert not selection["decision_label"].isna().any()
    card = json.loads(Path("reports/final_candidate_dossiers/cards/fcand_001.json").read_text(encoding="utf-8"))
    validate_candidate_card(card)
