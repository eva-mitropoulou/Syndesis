from __future__ import annotations

RESCORING_TASK_COLUMNS = [
    "rescoring_task_id", "pose_id", "docking_task_id", "ligand_id", "target_receptor_id",
    "task_type", "docking_engine", "original_pose_rank", "original_docking_score",
    "pose_file", "receptor_file", "receptor_state", "native_like_label_stage3", "rmsd_symmetry_corrected",
    "sanity_status", "rescoring_engine", "rescoring_mode", "model_name", "model_version",
    "use_gpu", "batch_id", "task_status", "skip_reason",
]

GNINA_SCORE_COLUMNS = [
    "pose_id", "rescoring_task_id", "docking_task_id", "ligand_id", "target_receptor_id",
    "task_type", "docking_engine", "original_pose_rank", "original_docking_score",
    "rmsd_symmetry_corrected", "sanity_status", "stage3_pose_label", "gnina_version",
    "gnina_model", "gnina_mode", "gnina_empirical_affinity", "cnnscore", "cnnaffinity",
    "cnn_vs", "cnnscore_rank_within_task", "cnnaffinity_rank_within_task",
    "gnina_affinity_rank_within_task", "rescoring_status", "rescoring_warnings_json",
]

RAW_GNINA_COLUMNS = [
    "rescoring_task_id", "pose_id", "gnina_version", "gnina_model", "gnina_mode",
    "command_line", "runtime_seconds", "exit_code", "status", "stdout_log", "stderr_log",
    "error_message",
]

EMPIRICAL_SCORE_COLUMNS = [
    "pose_id", "original_docking_score", "vina_rescore", "vinardo_rescore",
    "gnina_empirical_affinity", "empirical_rescoring_status", "warnings_json",
]

POSE_SCORE_COLUMNS = [
    "pose_id", "docking_task_id", "ligand_id", "target_receptor_id", "task_type",
    "receptor_state", "docking_engine", "original_pose_rank", "original_docking_score",
    "vina_rescore", "vinardo_rescore", "gnina_empirical_affinity", "cnnscore",
    "cnnaffinity", "cnn_vs", "rmsd_symmetry_corrected", "sanity_status",
    "stage3_pose_label", "strict_native_like_flag", "relaxed_native_like_flag",
    "invalid_pose_flag", "sampled_not_ranked_flag", "ranking_failure_flag",
    "sampling_failure_flag", "interaction_recovery_status", "final_pose_label_status",
    "rescoring_status", "warnings_json",
]

TASK_METRIC_COLUMNS = [
    "docking_task_id", "task_type", "ligand_id", "target_receptor_id", "docking_engine",
    "num_poses_scored", "top1_rmsd_original_score", "top1_rmsd_cnnscore",
    "top1_rmsd_cnnaffinity", "top1_rmsd_gnina_empirical", "best_rmsd_top5_original_score",
    "best_rmsd_top5_cnnscore", "best_rmsd_top5_cnnaffinity", "best_rmsd_top10_original_score",
    "best_rmsd_top10_cnnscore", "best_rmsd_top10_cnnaffinity",
    "top1_success_original_score_strict", "top1_success_cnnscore_strict",
    "top1_success_cnnaffinity_strict", "top1_success_original_score_relaxed",
    "top1_success_cnnscore_relaxed", "top1_success_cnnaffinity_relaxed",
    "rank_of_first_strict_native_like_by_original_score",
    "rank_of_first_strict_native_like_by_cnnscore",
    "rank_of_first_strict_native_like_by_cnnaffinity", "score_native_like_auc_original",
    "score_native_like_auc_cnnscore", "score_native_like_auc_cnnaffinity",
    "spearman_score_vs_rmsd_original", "spearman_score_vs_rmsd_cnnscore",
    "spearman_score_vs_rmsd_cnnaffinity", "dominant_best_scorer", "warnings_json",
]

SCORER_SUMMARY_COLUMNS = [
    "scorer_name", "num_tasks_evaluated", "redocking_top1_success_strict",
    "redocking_top1_success_relaxed", "same_state_crossdock_top1_success_strict",
    "same_state_crossdock_top1_success_relaxed", "other_state_crossdock_top1_success_strict",
    "other_state_crossdock_top1_success_relaxed", "mean_rank_of_first_strict_native_like",
    "median_rank_of_first_strict_native_like", "score_native_like_auc_mean",
    "score_native_like_auc_median", "spearman_score_vs_rmsd_median",
    "invalid_pose_preference_rate", "missing_score_rate", "runtime_seconds_total",
    "runtime_seconds_per_pose_median",
]

FAILURE_COLUMNS = [
    "pose_id", "docking_task_id", "ligand_id", "target_receptor_id", "task_type",
    "failure_type", "failure_description", "original_docking_score", "cnnscore",
    "cnnaffinity", "rmsd_symmetry_corrected", "sanity_status", "stage3_pose_label",
    "suggested_followup",
]

STAGE6_FEATURE_COLUMNS = [
    "pose_id", "docking_task_id", "ligand_id", "target_receptor_id", "task_type",
    "receptor_state", "docking_engine", "original_pose_rank", "original_docking_score",
    "vina_rescore", "vinardo_rescore", "gnina_empirical_affinity", "cnnscore",
    "cnnaffinity", "cnn_vs", "stage3_pose_label", "strict_native_like_flag",
    "relaxed_native_like_flag", "invalid_pose_flag", "rmsd_symmetry_corrected",
    "sanity_status",
]
