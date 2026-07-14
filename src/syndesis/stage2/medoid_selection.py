from __future__ import annotations

from typing import Any

import pandas as pd


def select_ensemble(features: pd.DataFrame, clusters: pd.DataFrame, config: dict[str, Any]) -> tuple[pd.DataFrame, pd.DataFrame]:
    merged = features.merge(clusters, on=["receptor_id", "state_stratum"], how="left")
    force_exclude = set(config["selection"].get("force_exclude", []))
    force_include = set(config["selection"].get("force_include_if_available", []))
    max_size = int(
        config["selection"].get(
            "preferred_ensemble_size_max",
            config["selection"].get("target_ensemble_size_max", 12),
        )
    )
    min_size = int(config["selection"].get("target_ensemble_size_min", 4))
    max_per_state = int(config["selection"].get("max_receptors_per_state_stratum", max_size))
    state_counts: dict[str, int] = {}

    candidates = merged[~merged["pdb_id"].isin(force_exclude)].copy()
    candidates["selection_rank"] = (
        candidates["cluster_medoid_flag"].fillna(False).astype(int) * 1000
        + candidates["receptor_preselection_score"].fillna(0)
        - candidates["distance_to_cluster_medoid"].fillna(0)
    )
    selected_ids: list[str] = []
    reasons: dict[str, str] = {}

    for pdb_id in force_include:
        control_rows = candidates[candidates["pdb_id"] == pdb_id].sort_values("selection_rank", ascending=False)
        if not control_rows.empty:
            rid = str(control_rows.iloc[0]["receptor_id"])
            selected_ids.append(rid)
            state = str(control_rows.iloc[0].get("state_stratum") or "unknown_state")
            state_counts[state] = state_counts.get(state, 0) + 1
            reasons[rid] = f"Reference control {pdb_id} passed Stage 2 filters."

    medoids = candidates[candidates["cluster_medoid_flag"] == True].sort_values("selection_rank", ascending=False)
    for _, row in medoids.iterrows():
        rid = str(row["receptor_id"])
        state = str(row.get("state_stratum") or "unknown_state")
        if state_counts.get(state, 0) >= max_per_state:
            continue
        if rid not in selected_ids:
            selected_ids.append(rid)
            state_counts[state] = state_counts.get(state, 0) + 1
            reasons[rid] = "Cluster medoid selected for pocket/state diversity."
        if len(selected_ids) >= max_size:
            break

    if len(selected_ids) < min_size:
        remaining = candidates[~candidates["receptor_id"].isin(selected_ids)].sort_values("selection_rank", ascending=False)
        for _, row in remaining.iterrows():
            rid = str(row["receptor_id"])
            state = str(row.get("state_stratum") or "unknown_state")
            selected_ids.append(rid)
            state_counts[state] = state_counts.get(state, 0) + 1
            reasons[rid] = "Added to satisfy minimum ensemble size."
            if len(selected_ids) >= min_size:
                break

    selected = merged[merged["receptor_id"].isin(selected_ids)].copy()
    selected["selected_flag"] = True
    selected["selected_role"] = selected["receptor_id"].map(
        lambda rid: "reference_control" if selected[selected["receptor_id"] == rid]["pdb_id"].iloc[0] in force_include else "cluster_medoid"
    )
    selected["selected_reason"] = selected["receptor_id"].map(reasons).fillna("Selected by medoid policy.")

    holdout = merged[~merged["receptor_id"].isin(selected_ids)].copy()
    holdout["selected_flag"] = False
    return selected, holdout
