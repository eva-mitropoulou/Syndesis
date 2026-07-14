from __future__ import annotations

from pathlib import Path

import pytest

from syndesis.stage12.candidate_dossiers import run_stage12_all
from syndesis.stage12.nonclaim_generator import FORBIDDEN_PHRASES, REQUIRED_NON_CLAIMS


CONFIG = Path("configs/stage12_candidate_dossiers.yaml")
pytestmark = pytest.mark.integration


def test_reports_contain_required_nonclaims_and_no_forbidden_phrases() -> None:
    run_stage12_all(CONFIG)
    report_paths = [Path("reports/12_final_candidate_summary.html"), *Path("reports/final_candidate_dossiers/html").glob("*.html")]
    for path in report_paths:
        text = path.read_text(encoding="utf-8")
        for claim in REQUIRED_NON_CLAIMS:
            assert claim in text
        lowered = text.lower()
        for phrase in FORBIDDEN_PHRASES:
            assert phrase not in lowered
