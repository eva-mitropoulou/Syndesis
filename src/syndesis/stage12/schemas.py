from __future__ import annotations

FINAL_CANDIDATE_COLUMNS = [
    "final_candidate_id", "molecule_id", "analog_id_if_available", "source", "subsource", "screening_role",
    "standard_smiles", "scaffold_id", "novelty_bucket", "closest_known_molecule_id",
    "tanimoto_to_closest_known", "best_pose_id", "best_receptor_id", "best_receptor_state",
    "final_candidate_score", "best_pose_confidence", "best_calibrated_confidence", "best_gnina_cnnscore",
    "best_gnina_cnnaffinity", "best_docking_score", "best_key_interaction_recall_consensus",
    "best_ifp_tanimoto_to_consensus", "md_stability_label_if_available",
    "md_key_interaction_persistence_if_available", "medchem_risk_score", "risk_flags_json",
    "decision_label", "selected_for_detailed_dossier", "selected_for_summary_table", "selection_reason",
    "nonclaim_statement", "warnings_json",
]

FINAL_RANKED_COLUMNS = [
    "final_rank_global", "final_rank_within_source", "final_rank_within_novelty_bucket",
    "final_candidate_id", "molecule_id", "source", "novelty_bucket", "standard_smiles",
    "closest_known_molecule_id", "tanimoto_to_closest_known", "best_receptor_state",
    "final_candidate_score", "pose_confidence", "calibrated_confidence", "gnina_cnnscore",
    "cnnaffinity", "key_interaction_recall_consensus", "ifp_tanimoto_to_consensus",
    "md_stability_label_if_available", "medchem_risk_score", "decision_label",
    "selected_for_dossier", "recommended_next_action", "nonclaim_statement",
]

ALLOWED_DECISION_LABELS = {
    "strong_computational_candidate", "binding_mode_preserved_analog", "high_confidence_close_analog",
    "high_confidence_medium_similarity", "scaffold_novel_but_risky", "generated_candidate_promising",
    "manual_candidate_promising", "known_control_recovered", "score_hacking_rejected",
    "md_unstable_rejected", "low_confidence_rejected", "pending_manual_review",
}


def validate_final_candidates(frame) -> None:
    missing = [column for column in FINAL_CANDIDATE_COLUMNS if column not in frame.columns]
    if missing:
        raise ValueError(f"Missing Stage 12 final candidate columns: {missing}")
    bad_labels = sorted(set(frame["decision_label"]) - ALLOWED_DECISION_LABELS)
    if bad_labels:
        raise ValueError(f"Unexpected decision labels: {bad_labels}")
