from __future__ import annotations

from pathlib import Path

import pytest

from egfr_dockingforge.stage12.report_stage12 import REQUIRED_FIGURES, render_final_figures


CONFIG = Path("configs/stage12_candidate_dossiers.yaml")
pytestmark = pytest.mark.integration


def test_required_stage12_figures_are_generated_and_referenced() -> None:
    render_final_figures(CONFIG)
    for filename in REQUIRED_FIGURES.values():
        path = Path("figures/stage12") / filename
        assert path.exists()
        assert path.stat().st_size > 100
