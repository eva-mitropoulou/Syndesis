from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pandas as pd

from egfr_dockingforge.common.io import write_table


def select_seed_scaffolds(inputs: dict[str, pd.DataFrame], config: dict[str, Any], paths: dict[str, Path]) -> pd.DataFrame:
    ranked = inputs["ranked_candidates"].copy()
    agg = inputs["candidate_aggregate_scores"].copy()
    if "candidate_rank" not in ranked.columns:
        ranked = ranked.sort_values("final_candidate_score", ascending=False).reset_index(drop=True)
        ranked["candidate_rank"] = ranked.index + 1
    # Pull best_gnina_cnnaffinity from the aggregate scores: the ranked table
    # only carries cnnscore, but parent ligand efficiency must be computed from
    # cnnaffinity to match the analog side (see analog_acceptance).
    agg_cols = ["molecule_id", "best_screening_pose_id", "best_target_receptor_id", "best_receptor_state"]
    if "best_gnina_cnnaffinity" in agg.columns:
        agg_cols.append("best_gnina_cnnaffinity")
    df = ranked.merge(
        agg[agg_cols],
        on="molecule_id",
        how="left",
    )
    filt = df[
        (df["best_pose_confidence"] >= float(config["seed_selection"]["min_pose_confidence"]))
        & (df["best_key_interaction_recall_consensus"] >= float(config["seed_selection"]["min_key_interaction_recall"]))
    ].copy()
    if filt.empty:
        filt = df.head(int(config["seed_selection"]["min_seeds"])).copy()
        reason = "top_ranked_stage8_control_seed_relaxed_due_to_no_prospective_candidates"
    else:
        reason = "passed_stage8_pose_confidence_and_interaction_filters"
    filt = filt.head(int(config["seed_selection"]["max_seeds"]))
    rows = []
    for i, row in enumerate(filt.to_dict("records"), start=1):
        rows.append(
            {
                "seed_id": f"seed_{i:03d}",
                "molecule_id": row["molecule_id"],
                "source": row.get("source", ""),
                "standard_smiles": row["standard_smiles"],
                "scaffold_id": row.get("scaffold_id", ""),
                "parent_candidate_rank": int(row.get("candidate_rank", i)),
                "best_pose_id": row.get("best_screening_pose_id", row.get("best_pose_id", "")),
                "best_receptor_id": row.get("best_target_receptor_id", row.get("best_receptor_id", "")),
                "best_receptor_state": row.get("best_receptor_state", ""),
                "best_pose_confidence": float(row.get("best_pose_confidence", 0.0)),
                "best_gnina_cnnscore": float(row.get("best_gnina_cnnscore", 0.0)),
                "best_gnina_cnnaffinity": float(row.get("best_gnina_cnnaffinity", 0.0)),
                "best_key_interaction_recall_consensus": float(row.get("best_key_interaction_recall_consensus", 0.0)),
                "novelty_bucket": row.get("novelty_bucket", "known_control"),
                "medchem_flags_json": row.get("medchem_flags_json", "[]"),
                "selected_for_analog_loop": True,
                "selection_reason": reason,
                "warnings_json": json.dumps(["stage8_input_contains_known_controls_only"] if row.get("screening_role") == "known_activity_reference" else []),
            }
        )
    out = pd.DataFrame(rows)
    write_table(paths["processed"] / "analog_seed_scaffolds.parquet", out)
    write_table(paths["processed"] / "analog_seed_scaffolds.csv", out)
    return out
