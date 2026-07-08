from __future__ import annotations

from typing import Any

import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import squareform

from egfr_dockingforge.stage2.pocket_alignment import receptor_pair_rmsd
from egfr_dockingforge.stage2.schemas import CLUSTER_COLUMNS


def distance_matrix(features: pd.DataFrame, config: dict[str, Any]) -> pd.DataFrame:
    residues = [int(x) for x in config["pocket"]["alignment_residues_uniprot"]]
    ids = list(features["receptor_id"])
    rows = []
    for i, id_a in enumerate(ids):
        path_a = features.iloc[i]["receptor_file_path"]
        for j, id_b in enumerate(ids):
            path_b = features.iloc[j]["receptor_file_path"]
            if i == j:
                distance = 0.0
            elif j < i:
                continue
            else:
                value, _, _ = receptor_pair_rmsd(path_a, path_b, residues)
                distance = float(value) if value is not None else 999.0
            rows.append({"receptor_id_a": id_a, "receptor_id_b": id_b, "pocket_ca_rmsd": distance})
            if j > i:
                rows.append({"receptor_id_a": id_b, "receptor_id_b": id_a, "pocket_ca_rmsd": distance})
    return pd.DataFrame(rows)


def matrix_for_ids(distance_table: pd.DataFrame, ids: list[str]) -> np.ndarray:
    lookup = {
        (row.receptor_id_a, row.receptor_id_b): float(row.pocket_ca_rmsd)
        for row in distance_table.itertuples(index=False)
    }
    matrix = np.zeros((len(ids), len(ids)), dtype=float)
    for i, a in enumerate(ids):
        for j, b in enumerate(ids):
            matrix[i, j] = lookup.get((a, b), 999.0)
    return matrix


def medoid_for_cluster(ids: list[str], matrix: np.ndarray, members: list[int]) -> tuple[str, dict[str, float]]:
    sub = matrix[np.ix_(members, members)]
    sums = sub.sum(axis=1)
    best_local = int(np.argmin(sums))
    medoid_index = members[best_local]
    medoid = ids[medoid_index]
    distances = {ids[members[i]]: float(sub[i, best_local]) for i in range(len(members))}
    return medoid, distances


def cluster_receptors(features: pd.DataFrame, distance_table: pd.DataFrame, config: dict[str, Any]) -> pd.DataFrame:
    threshold = float(config["clustering"].get("state_stratified_pocket_rmsd_threshold", 1.25))
    rows = []
    for state, group in features.groupby("state_stratum", dropna=False):
        ids = list(group["receptor_id"])
        if len(ids) == 1:
            rows.append(
                {
                    "receptor_id": ids[0],
                    "state_stratum": state or "unknown_state",
                    "cluster_id": f"{state}_1",
                    "cluster_size": 1,
                    "cluster_medoid_flag": True,
                    "distance_to_cluster_medoid": 0.0,
                    "nearest_neighbor_receptor_id": None,
                    "nearest_neighbor_distance": None,
                    "cluster_quality_summary": "singleton",
                    "cluster_warning": "",
                }
            )
            continue
        matrix = matrix_for_ids(distance_table, ids)
        condensed = squareform(matrix, checks=False)
        labels = fcluster(linkage(condensed, method="average"), t=threshold, criterion="distance")
        for label in sorted(set(labels)):
            member_indices = [i for i, value in enumerate(labels) if value == label]
            medoid, distances = medoid_for_cluster(ids, matrix, member_indices)
            for idx in member_indices:
                receptor_id = ids[idx]
                nearest = None
                nearest_distance = None
                if len(ids) > 1:
                    row_dist = matrix[idx].copy()
                    row_dist[idx] = np.inf
                    nearest_idx = int(np.argmin(row_dist))
                    nearest = ids[nearest_idx]
                    nearest_distance = float(row_dist[nearest_idx])
                rows.append(
                    {
                        "receptor_id": receptor_id,
                        "state_stratum": state or "unknown_state",
                        "cluster_id": f"{state}_{label}",
                        "cluster_size": len(member_indices),
                        "cluster_medoid_flag": receptor_id == medoid,
                        "distance_to_cluster_medoid": distances[receptor_id],
                        "nearest_neighbor_receptor_id": nearest,
                        "nearest_neighbor_distance": nearest_distance,
                        "cluster_quality_summary": f"{len(member_indices)} receptors",
                        "cluster_warning": "",
                    }
                )
    return pd.DataFrame(rows, columns=CLUSTER_COLUMNS)

