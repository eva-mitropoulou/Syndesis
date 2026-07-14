from __future__ import annotations

RECEPTOR_PREP_COLUMNS = [
    "receptor_id", "complex_id", "pdb_id", "auth_asym_id", "receptor_state",
    "input_receptor_file", "prepared_receptor_file", "docking_format_file", "engine",
    "hydrogens_added_flag", "protonation_tool", "retained_waters_flag", "retained_water_count",
    "preparation_warnings_json", "preparation_status",
]

LIGAND_PREP_COLUMNS = [
    "ligand_id", "complex_id", "pdb_id_native", "ligand_comp_id", "ligand_instance_id",
    "native_ligand_file", "immutable_reference_pose_file", "standard_smiles", "prepared_smiles",
    "protonation_state_id", "tautomer_state_id", "conformer_id", "prepared_ligand_file",
    "docking_format_file", "engine", "ligand_prep_tool", "charge_model", "atom_mapping_status",
    "native_to_prepared_atom_map_json", "preparation_warnings_json", "preparation_status",
]

TASK_COLUMNS = [
    "docking_task_id", "task_type", "ligand_id", "ligand_native_complex_id", "native_receptor_id",
    "target_receptor_id", "target_receptor_state", "native_receptor_state", "state_match_flag",
    "ligand_prepared_file", "receptor_prepared_file", "docking_box_center_x", "docking_box_center_y",
    "docking_box_center_z", "docking_box_size_x", "docking_box_size_y", "docking_box_size_z",
    "docking_engine", "engine_version", "exhaustiveness", "num_modes", "seed", "replicate_id",
    "expected_native_pose_reference_file", "reference_pose_transform_matrix_path", "task_status",
    "skip_reason",
]

TRANSFORM_COLUMNS = [
    "ligand_id", "native_receptor_id", "target_receptor_id", "transform_matrix_file",
    "pocket_alignment_rmsd", "alignment_residue_count", "transformed_reference_pose_file",
    "transform_status", "warnings_json",
]

DOCKING_RUN_COLUMNS = [
    "docking_run_id", "docking_task_id", "docking_engine", "engine_version", "command_line",
    "container_image", "start_time", "end_time", "runtime_seconds", "exit_code", "status",
    "error_message", "output_pose_file", "output_log_file", "config_hash",
]

DOCKED_POSE_COLUMNS = [
    "pose_id", "docking_run_id", "docking_task_id", "ligand_id", "target_receptor_id",
    "docking_engine", "protonation_state_id", "tautomer_state_id", "conformer_id", "seed",
    "replicate_id", "pose_rank", "docking_score", "pose_file", "raw_pose_file", "parse_status",
    "parse_warnings_json",
]

RMSD_COLUMNS = [
    "pose_id", "docking_task_id", "ligand_id", "target_receptor_id", "task_type", "docking_engine",
    "pose_rank", "docking_score", "rmsd_heavy_atom", "rmsd_symmetry_corrected", "rmsd_method",
    "atom_mapping_status", "atom_mapping_warning", "reference_pose_file", "transformed_reference_pose_flag",
    "pocket_alignment_rmsd_for_reference", "rmsd_status", "rmsd_warnings_json",
]

SANITY_COLUMNS = [
    "pose_id", "sanity_status", "severe_clash_flag", "ligand_geometry_flag", "chirality_issue_flag",
    "atom_loss_flag", "outside_pocket_flag", "protein_ligand_clash_score", "ligand_strain_proxy",
    "posebusters_available_flag", "posebusters_pass_flag", "sanity_warnings_json",
]

LABEL_COLUMNS = [
    "pose_id", "docking_task_id", "ligand_id", "target_receptor_id", "task_type", "docking_engine",
    "pose_rank", "docking_score", "rmsd_symmetry_corrected", "sanity_status",
    "strict_native_like_flag", "relaxed_native_like_flag", "invalid_pose_flag",
    "sampled_not_ranked_flag", "ranking_failure_flag", "sampling_failure_flag", "stage3_pose_label",
    "interaction_recovery_status", "final_pose_label_status", "label_warnings_json",
]

TASK_METRIC_COLUMNS = [
    "docking_task_id", "task_type", "ligand_id", "native_receptor_id", "target_receptor_id",
    "state_match_flag", "docking_engine", "num_poses_generated", "top1_rmsd", "top1_score",
    "top5_best_rmsd", "top10_best_rmsd", "top20_best_rmsd", "best_rmsd_any_pose",
    "rank_of_first_strict_native_like", "rank_of_first_relaxed_native_like",
    "sampling_success_strict_top20", "sampling_success_relaxed_top20",
    "ranking_success_strict_top1", "ranking_success_relaxed_top1", "physical_sanity_pass_rate",
    "failed_run_flag", "failure_category", "warnings_json",
]

RECEPTOR_VALIDATION_COLUMNS = [
    "receptor_id", "pdb_id", "receptor_state", "num_redocking_tasks", "num_crossdocking_tasks",
    "redocking_top1_success_rate_strict", "redocking_top20_success_rate_strict",
    "same_state_crossdock_top1_success_rate_strict", "same_state_crossdock_top20_success_rate_strict",
    "other_state_crossdock_top20_success_rate_relaxed", "physical_sanity_pass_rate",
    "dominant_failure_mode", "keep_for_stage4_flag", "prune_recommendation_flag",
    "recommendation_reason", "warnings_json",
]
