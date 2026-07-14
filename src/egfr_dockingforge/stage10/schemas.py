from __future__ import annotations

TABLE_SCHEMAS = {
    "ablation_strategy_manifest": [
        "strategy_id","strategy_name","strategy_family","uses_rdkit_validation","uses_docking_score","uses_gnina","uses_prolif_constraint","uses_pose_confidence","uses_md_filter","uses_medchem_filter","generation_budget","screening_budget","seed_set_id","num_seed_scaffolds","num_iterations","enabled_flag","md_status","notes",
    ],
    "strategy_budget_audit": [
        "strategy_id","seed_id","num_raw_proposals","num_valid_unique_analogs","num_screened_analogs","num_docking_tasks","num_gnina_tasks","num_prolif_tasks","num_pose_confidence_predictions","num_md_tasks_if_available","walltime_seconds","gpu_seconds","cpu_seconds","budget_normalized_flag","budget_violation_reason",
    ],
    "analog_benchmark_master": [
        "analog_id","strategy_id","strategy_name","seed_id","parent_molecule_id","parent_smiles","analog_smiles","standard_smiles","inchi_key","scaffold_id","valid_molecule_flag","unique_flag","hard_scope_pass","covalent_warhead_flag","reactive_flag","pains_flag","property_window_pass","medchem_risk_score","parent_tanimoto","closest_known_egfr_ligand","novelty_bucket","best_docking_score","best_gnina_cnnscore","best_gnina_cnnaffinity","best_pose_confidence","best_key_interaction_recall_consensus","best_ifp_tanimoto_to_consensus","binding_mode_preserved_flag","ligand_efficiency","parent_candidate_score","analog_candidate_score","delta_candidate_score","delta_pose_confidence","delta_gnina_cnnscore","delta_key_interaction_recall","delta_ligand_efficiency","md_stability_label_if_available","md_key_interaction_persistence_if_available","accepted_pre_md_flag","accepted_post_md_flag","rejection_reason","score_hacking_flag","warnings_json",
    ],
    "score_hacking_cases": ["analog_id","strategy_id","seed_id","score_hacking_type","improved_metric","worsened_metric","parent_value","analog_value","severity","evidence_json","warnings_json"],
    "strategy_metrics": [
        "strategy_id","strategy_name","num_seeds","num_raw_proposals","num_valid_molecules","num_unique_molecules","num_screened","num_pre_md_accepted","num_post_md_accepted_if_available","validity_rate","uniqueness_rate","novelty_rate","diversity_score","accepted_analog_rate_pre_md","accepted_analog_rate_post_md_if_available","score_hacking_rate","binding_mode_break_rate","bad_chemistry_rejection_rate","medchem_risk_rejection_rate","mean_delta_candidate_score","median_delta_candidate_score","mean_delta_pose_confidence","median_delta_pose_confidence","mean_delta_key_interaction_recall","median_delta_key_interaction_recall","mean_delta_ligand_efficiency","median_delta_ligand_efficiency","mean_runtime_per_valid_analog","mean_runtime_per_accepted_analog","gpu_seconds_per_accepted_analog","md_status","notes",
    ],
    "seed_strategy_metrics": [
        "seed_id","strategy_id","strategy_name","num_raw_proposals","num_valid_molecules","num_unique_molecules","num_screened","num_accepted_pre_md","num_accepted_post_md_if_available","accepted_rate_pre_md","accepted_rate_post_md_if_available","score_hacking_rate","best_delta_candidate_score","best_delta_pose_confidence","best_delta_key_interaction_recall","best_delta_ligand_efficiency","best_accepted_analog_id","dominant_rejection_reason","warnings_json",
    ],
    "statistical_comparisons": ["comparison_id","metric_name","method_a","method_b","mean_a","mean_b","median_a","median_b","delta_mean","effect_size","ci_low","ci_high","p_value","p_value_corrected","test_name","num_seed_pairs","interpretation","warnings_json"],
    "ablation_summary": ["ablation_id","base_strategy","removed_component","added_component","accepted_rate_change","score_hacking_rate_change","binding_mode_break_rate_change","medchem_risk_change","runtime_change","conclusion","evidence_strength","warnings_json"],
    "diversity_novelty_metrics": ["strategy_id","seed_id","internal_diversity","unique_scaffold_count","unique_scaffold_rate","mean_parent_tanimoto","median_parent_tanimoto","known_duplicate_rate","close_analog_rate","scaffold_novel_rate","mode_collapse_flag","warnings_json"],
    "compute_cost_metrics": ["strategy_id","cpu_seconds_total","gpu_seconds_total","walltime_seconds_total","num_docking_tasks","num_gnina_tasks","num_prolif_tasks","num_md_tasks_if_available","accepted_analogs_per_gpu_hour","accepted_analogs_per_wall_hour","notes"],
}
