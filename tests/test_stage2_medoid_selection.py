from __future__ import annotations

import pandas as pd

from egfr_dockingforge.stage2.medoid_selection import select_ensemble


def test_medoid_selection_is_reproducible() -> None:
    features = pd.DataFrame(
        [
            {"receptor_id": "r1", "pdb_id": "1M17", "receptor_preselection_score": 90.0, "state_stratum": "active_like"},
            {"receptor_id": "r2", "pdb_id": "1XKK", "receptor_preselection_score": 95.0, "state_stratum": "inactive_like"},
        ]
    )
    clusters = pd.DataFrame(
        [
            {"receptor_id": "r1", "state_stratum": "active_like", "cluster_id": "active_like_1", "cluster_medoid_flag": True, "distance_to_cluster_medoid": 0.0},
            {"receptor_id": "r2", "state_stratum": "inactive_like", "cluster_id": "inactive_like_1", "cluster_medoid_flag": True, "distance_to_cluster_medoid": 0.0},
        ]
    )
    config = {
        "selection": {
            "force_exclude": ["4ZAU"],
            "force_include_if_available": ["1M17", "1XKK"],
            "target_ensemble_size_max": 12,
            "target_ensemble_size_min": 2,
        }
    }
    selected_a, _ = select_ensemble(features, clusters, config)
    selected_b, _ = select_ensemble(features, clusters, config)
    assert list(selected_a["receptor_id"]) == list(selected_b["receptor_id"])
    assert set(selected_a["receptor_id"]) == {"r1", "r2"}

