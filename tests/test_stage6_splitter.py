from __future__ import annotations

import pandas as pd

from egfr_dockingforge.stage6.splitter import build_ranking_groups


def test_ranking_groups_exclude_no_variation_groups(tmp_path) -> None:
    features = pd.DataFrame(
        {
            "pose_id": ["p1", "p2", "p3", "p4"],
            "docking_task_id": ["g1", "g1", "g2", "g2"],
            "ligand_id": ["l1", "l1", "l2", "l2"],
            "target_receptor_id": ["r1", "r1", "r2", "r2"],
        }
    )
    labels = pd.DataFrame({"pose_id": ["p1", "p2", "p3", "p4"], "rank_relevance_label": [0, 1, 0, 0]})
    groups = build_ranking_groups(features, labels, {"processed": tmp_path})
    assert groups.loc[groups["group_id"].eq("g1"), "group_usable_for_ranking"].iloc[0]
    assert not groups.loc[groups["group_id"].eq("g2"), "group_usable_for_ranking"].iloc[0]
