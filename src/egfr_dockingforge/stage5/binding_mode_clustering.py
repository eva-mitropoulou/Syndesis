from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import pdist, squareform

from egfr_dockingforge.common.io import write_table
from egfr_dockingforge.stage5.schemas import BINDING_MODE_CLUSTER_COLUMNS


def _bits(value: object) -> set[str]:
    try:
        return set(json.loads(value if isinstance(value, str) and value else "[]"))
    except json.JSONDecodeError:
        return set()


def _jaccard_distance(a: set[str], b: set[str]) -> float:
    union = a | b
    if not union:
        return 1.0
    return 1.0 - len(a & b) / len(union)


def cluster_binding_modes(
    native_fps: pd.DataFrame,
    docked_fps: pd.DataFrame,
    inputs: dict[str, pd.DataFrame],
    config: dict[str, Any],
    paths: dict[str, Path],
) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for row in native_fps.itertuples(index=False):
        rows.append(
            {
                "entity_id": row.complex_id,
                "entity_type": "native_complex",
                "pose_id": None,
                "complex_id": row.complex_id,
                "ligand_id": row.ligand_id,
                "receptor_id": row.receptor_id,
                "receptor_state": row.receptor_state,
                "fingerprint_bitstring": row.fingerprint_bitstring,
                "_bits": _bits(row.fingerprint_sparse_json),
            }
        )
    receptor_state = inputs["pose_scores"].set_index("pose_id")["receptor_state"].to_dict()
    for row in docked_fps.itertuples(index=False):
        rows.append(
            {
                "entity_id": row.pose_id,
                "entity_type": "docked_pose",
                "pose_id": row.pose_id,
                "complex_id": None,
                "ligand_id": row.ligand_id,
                "receptor_id": row.target_receptor_id,
                "receptor_state": receptor_state.get(row.pose_id),
                "fingerprint_bitstring": row.fingerprint_bitstring,
                "_bits": _bits(row.fingerprint_sparse_json),
            }
        )
    if not rows:
        frame = pd.DataFrame(columns=BINDING_MODE_CLUSTER_COLUMNS)
        write_table(paths["processed"] / "binding_mode_clusters.parquet", frame)
        write_table(paths["processed"] / "binding_mode_clusters.csv", frame)
        return frame
    bits = [row["_bits"] for row in rows]
    if len(rows) == 1:
        labels = np.array([1])
        distances = np.array([[0.0]])
    else:
        condensed = pdist(np.array(range(len(bits))).reshape(-1, 1), metric=lambda u, v: _jaccard_distance(bits[int(u[0])], bits[int(v[0])]))
        labels = fcluster(linkage(condensed, method="average"), t=float(config.get("clustering", {}).get("distance_threshold", 0.65)), criterion="distance")
        distances = squareform(condensed)
    cluster_sizes = pd.Series(labels).value_counts().to_dict()
    medoids: dict[int, int] = {}
    for cluster_id in sorted(set(labels)):
        members = [idx for idx, label in enumerate(labels) if label == cluster_id]
        sub = distances[np.ix_(members, members)]
        medoids[int(cluster_id)] = members[int(np.argmin(sub.mean(axis=1)))]
    out_rows: list[dict[str, Any]] = []
    for idx, row in enumerate(rows):
        cluster_id = int(labels[idx])
        medoid = medoids[cluster_id]
        dominant = sorted(set().union(*[bits[i] for i, label in enumerate(labels) if int(label) == cluster_id]))
        out_rows.append(
            {
                **{k: v for k, v in row.items() if k != "_bits"},
                "cluster_id": f"ifp_cluster_{cluster_id}",
                "cluster_label": "native_enriched" if any(rows[i]["entity_type"] == "native_complex" for i, label in enumerate(labels) if int(label) == cluster_id) else "docked_only",
                "cluster_size": int(cluster_sizes[cluster_id]),
                "cluster_medoid_flag": idx == medoid,
                "distance_to_medoid": float(distances[idx, medoid]),
                "dominant_key_interactions_json": json.dumps(dominant[:25]),
                "cluster_interpretation": "interaction-fingerprint cluster",
                "warnings_json": json.dumps([]),
            }
        )
    frame = pd.DataFrame(out_rows, columns=BINDING_MODE_CLUSTER_COLUMNS)
    write_table(paths["processed"] / "binding_mode_clusters.parquet", frame)
    write_table(paths["processed"] / "binding_mode_clusters.csv", frame)
    return frame
