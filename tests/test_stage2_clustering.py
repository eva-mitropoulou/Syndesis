from __future__ import annotations

import pandas as pd

from syndesis.stage2.pocket_clustering import cluster_receptors


def test_singleton_cluster_has_medoid() -> None:
    features = pd.DataFrame(
        [{"receptor_id": "r1", "state_stratum": "active_like", "receptor_file_path": "unused.pdb"}]
    )
    distances = pd.DataFrame([{"receptor_id_a": "r1", "receptor_id_b": "r1", "pocket_ca_rmsd": 0.0}])
    config = {"clustering": {"state_stratified_pocket_rmsd_threshold": 1.0}}
    clusters = cluster_receptors(features, distances, config)
    assert clusters.loc[0, "cluster_medoid_flag"] is True or clusters.loc[0, "cluster_medoid_flag"] == True
    assert clusters.loc[0, "cluster_id"] == "active_like_1"

