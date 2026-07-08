from __future__ import annotations

import pandas as pd

from egfr_dockingforge.stage7.diversity_selection import select_subsets


def test_diversity_selection_keeps_known_controls_separate() -> None:
    master = pd.DataFrame([{"molecule_id": "m1", "source": "chembl_known_ligand", "subsource": "chembl", "screening_role": "known_activity_reference", "include_in_screening_library": False, "novelty_bucket": "known_duplicate", "scaffold_id": "s1"}])
    subsets = select_subsets(master, {"diversity": {"smoke_test_size": 10, "dev_screen_size": 10, "main_screen_size": 10}})
    assert "known_controls" in set(subsets["screening_subset"])
