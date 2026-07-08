from __future__ import annotations

import json
import pandas as pd


def rank_candidates(aggregate: pd.DataFrame, paths: dict) -> tuple[pd.DataFrame, pd.DataFrame]:
    ranked = aggregate.sort_values("final_candidate_score", ascending=False).copy()
    ranked["final_rank_global"] = range(1, len(ranked) + 1)
    ranked["final_rank_within_source"] = ranked.groupby("source")["final_candidate_score"].rank(ascending=False, method="first").astype(int)
    ranked["final_rank_within_novelty_bucket"] = ranked.groupby("novelty_bucket")["final_candidate_score"].rank(ascending=False, method="first").astype(int)
    ranked["candidate_bucket"] = ranked.apply(lambda r: "known_control_recovered" if r["screening_role"] in {"known_activity_reference", "native_pose_reference"} and r["candidate_decision_label"].startswith("accepted") else "low_confidence_rejected", axis=1)
    ranked["recommended_next_stage"] = ranked["candidate_bucket"].map({"known_control_recovered": "manual inspection"}).fillna("reject")
    ranked["selection_reason"] = ranked["candidate_decision_reason"]
    ranked["risk_summary"] = ranked["medchem_flags_json"]
    out_cols = ["final_rank_global", "final_rank_within_source", "final_rank_within_novelty_bucket", "molecule_id", "source", "screening_role", "novelty_bucket", "standard_smiles", "best_screening_pose_id", "best_target_receptor_id", "best_receptor_state", "final_candidate_score", "best_pose_confidence", "best_gnina_cnnscore", "best_key_interaction_recall_consensus", "best_ifp_tanimoto_to_consensus", "medchem_risk_score", "tanimoto_to_closest_known", "closest_known_molecule_id", "scaffold_id", "candidate_bucket", "recommended_next_stage", "selection_reason", "risk_summary", "warnings_json"]
    out = ranked[out_cols]
    out.to_parquet(paths["processed"] / "ranked_candidates.parquet", index=False)
    out.to_csv(paths["processed"] / "ranked_candidates.csv", index=False)
    controls = ranked[ranked["screening_role"].isin(["known_activity_reference", "native_pose_reference"])].copy()
    controls["known_activity_status"] = "known_reference"
    controls["endpoint_type_if_known"] = ""
    controls["median_p_activity_if_known"] = None
    controls["best_candidate_score"] = controls["final_candidate_score"]
    controls["rank_within_known_controls"] = controls["best_candidate_score"].rank(ascending=False, method="first").astype(int)
    controls["rank_within_full_screen"] = controls["final_rank_global"]
    controls["recovered_expected_binding_mode_flag"] = controls["candidate_bucket"].eq("known_control_recovered")
    controls["diagnostic_status"] = controls["recovered_expected_binding_mode_flag"].map({True: "pass", False: "warning"})
    diag = controls[["molecule_id", "known_activity_status", "endpoint_type_if_known", "median_p_activity_if_known", "source", "best_pose_confidence", "best_candidate_score", "rank_within_known_controls", "rank_within_full_screen", "recovered_expected_binding_mode_flag", "diagnostic_status", "warnings_json"]]
    diag.to_parquet(paths["processed"] / "known_control_diagnostics.parquet", index=False)
    diag.to_csv(paths["processed"] / "known_control_diagnostics.csv", index=False)
    pd.DataFrame(columns=["iteration_id", "candidate_pool_size", "num_docked_this_iteration", "surrogate_model_type", "acquisition_function", "selected_molecule_ids_json", "top_score_recovery_estimate", "stopping_reason", "warnings_json"]).to_parquet(paths["processed"] / "active_learning_iterations.parquet", index=False)
    pd.DataFrame(columns=["molecule_id", "tool", "score", "warnings_json"]).to_parquet(paths["processed"] / "finalist_optional_ai_scores.parquet", index=False)
    return out, diag
