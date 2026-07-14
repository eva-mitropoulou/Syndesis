from __future__ import annotations

IDENTIFIER_COLUMNS = [
    "pose_id",
    "docking_task_id",
    "ligand_id",
    "target_receptor_id",
    "native_receptor_id",
    "task_type",
    "docking_engine",
    "receptor_state",
    "ligand_source",
    "pose_rank",
]

LABEL_COLUMNS = [
    "pose_id",
    "docking_task_id",
    "final_pose_label",
    "stage3_pose_label",
    "rmsd_symmetry_corrected",
    "sanity_status",
    "interaction_recovery_label",
    "rank_relevance_label",
    "binary_confidence_label",
    "label_confidence",
    "label_reason",
]

LEAKAGE_AUDIT_COLUMNS = [
    "feature_name",
    "feature_group",
    "allowed_for_training",
    "allowed_for_deployment",
    "leakage_risk",
    "reason",
    "action",
]

RANKING_GROUP_COLUMNS = [
    "group_id",
    "docking_task_id",
    "ligand_id",
    "target_receptor_id",
    "num_poses",
    "num_positive_or_relevant_poses",
    "max_relevance_label",
    "group_usable_for_ranking",
    "exclusion_reason",
]

SPLIT_COLUMNS = [
    "pose_id",
    "group_id",
    "ligand_id",
    "scaffold_id",
    "target_receptor_id",
    "receptor_state",
    "split_name",
    "split_fold",
    "train_valid_test",
    "split_reason",
]

METRIC_COLUMNS = [
    "model_id",
    "model_family",
    "split_name",
    "split_fold",
    "metric_name",
    "metric_value",
    "metric_context",
    "num_groups",
    "num_poses",
    "notes",
]

FORBIDDEN_FEATURES = {
    "rmsd_symmetry_corrected",
    "rmsd_heavy_atom",
    "strict_native_like_flag",
    "relaxed_native_like_flag",
    "final_pose_label",
    "stage3_pose_label",
    "ifp_tanimoto_to_native",
    "key_interaction_recall_native",
    "key_interaction_precision_native",
    "key_interaction_f1_native",
    "native_reference_complex_id",
    "native_reference_available_flag",
    "transformed_reference_pose_flag",
    # --- Target leakage: the pose-confidence label is derived from the Stage 5
    # interaction-consistency assessment (final_pose_label) and native-like RMSD
    # flags. Any feature that measures a pose's agreement with the native /
    # consensus interaction reference therefore encodes (a transform of) the
    # label and must NOT enter the training matrix. The raw per-pose ProLIF
    # fingerprint bits (ifp_bit_*) are legitimate: they describe the pose itself,
    # not its agreement with the label.
    "ifp_tanimoto_to_consensus",
    "key_interaction_recall_consensus",
    "key_interaction_precision_consensus",
    "key_interaction_f1_consensus",
    "interaction_recovery_label",
    "hinge_interaction_recovered_flag",
    "catalytic_lys_glu_region_consistent_flag",
    "gatekeeper_region_consistent_flag",
    "dfg_region_consistent_flag",
    "binding_mode_reference_id",
    "missing_key_interactions_json",
    "extra_key_interactions_json",
    "recovered_key_interactions_json",
    # Native-cluster membership and redock-vs-crossdock design flags encode
    # native-likeness / the experimental design rather than a deployable signal.
    "dominant_binding_mode_cluster_compatibility",
    "cluster_label",
    "state_match_flag",
}
