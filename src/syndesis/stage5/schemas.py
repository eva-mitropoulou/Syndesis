from __future__ import annotations

INTERACTION_COMPLEX_COLUMNS = [
    "complex_analysis_id", "pose_id", "complex_id", "ligand_id", "receptor_id", "pdb_id",
    "task_type", "is_native_complex", "protein_file", "ligand_file", "complex_file",
    "ligand_source_coordinate_type", "hydrogens_present_flag",
    "ligand_bond_orders_available_flag", "protonation_state_id", "tautomer_state_id",
    "complex_build_status", "warnings_json",
]

RESIDUE_MAP_COLUMNS = [
    "receptor_id", "pdb_id", "auth_asym_id", "residue_name", "auth_seq_id",
    "uniprot_residue_number", "klifs_position", "residue_role", "mapping_source",
    "mapping_confidence", "warnings_json",
]

NATIVE_INTERACTIONS_LONG_COLUMNS = [
    "native_interaction_id", "complex_id", "receptor_id", "ligand_id", "pdb_id",
    "ligand_comp_id", "receptor_state", "mutation_flag", "mutation_list",
    "ligand_class_if_known", "residue_name", "auth_seq_id", "uniprot_residue_number",
    "klifs_position", "residue_role", "interaction_type", "present",
    "atom_indices_ligand_json", "atom_indices_protein_json", "distance", "angle",
    "interaction_confidence", "warnings_json",
]

NATIVE_FINGERPRINT_COLUMNS = [
    "complex_id", "receptor_id", "ligand_id", "pdb_id", "receptor_state",
    "ligand_class_if_known", "fingerprint_bitstring", "fingerprint_sparse_json",
    "num_interactions", "num_key_interactions", "interaction_engine",
    "interaction_engine_version", "interaction_config_hash", "warnings_json",
]

KEY_INTERACTION_COLUMNS = [
    "key_interaction_id", "receptor_state_scope", "binding_mode_scope", "residue_name",
    "uniprot_residue_number", "klifs_position", "residue_role", "interaction_type",
    "native_frequency", "native_count", "native_total", "selection_reason", "weight",
    "manual_override_flag", "evidence_complex_ids_json", "warnings_json",
]

DOCKED_INTERACTIONS_LONG_COLUMNS = [
    "pose_interaction_id", "pose_id", "docking_task_id", "ligand_id", "target_receptor_id",
    "native_receptor_id", "task_type", "docking_engine", "pose_rank", "docking_score",
    "gnina_cnnscore", "gnina_cnnaffinity", "rmsd_symmetry_corrected", "stage3_pose_label",
    "sanity_status", "receptor_state", "residue_name", "auth_seq_id",
    "uniprot_residue_number", "klifs_position", "residue_role", "interaction_type",
    "present", "key_interaction_flag", "atom_indices_ligand_json",
    "atom_indices_protein_json", "distance", "angle", "interaction_confidence",
    "warnings_json",
]

DOCKED_FINGERPRINT_COLUMNS = [
    "pose_id", "docking_task_id", "ligand_id", "target_receptor_id", "task_type",
    "docking_engine", "pose_rank", "fingerprint_bitstring", "fingerprint_sparse_json",
    "num_interactions", "num_key_interactions", "interaction_engine",
    "interaction_engine_version", "interaction_config_hash", "warnings_json",
]

INTERACTION_RECOVERY_COLUMNS = [
    "pose_id", "docking_task_id", "ligand_id", "target_receptor_id", "task_type",
    "docking_engine", "pose_rank", "rmsd_symmetry_corrected", "stage3_pose_label",
    "sanity_status", "native_reference_available_flag", "native_reference_complex_id",
    "binding_mode_reference_id", "ifp_tanimoto_to_native", "ifp_tanimoto_to_consensus",
    "key_interaction_recall_native", "key_interaction_precision_native",
    "key_interaction_f1_native", "key_interaction_recall_consensus",
    "key_interaction_precision_consensus", "key_interaction_f1_consensus",
    "missing_key_interactions_json", "extra_key_interactions_json",
    "recovered_key_interactions_json", "hinge_interaction_recovered_flag",
    "catalytic_lys_glu_region_consistent_flag", "gatekeeper_region_consistent_flag",
    "dfg_region_consistent_flag", "interaction_recovery_label", "warnings_json",
]

BINDING_MODE_CLUSTER_COLUMNS = [
    "entity_id", "entity_type", "pose_id", "complex_id", "ligand_id", "receptor_id",
    "receptor_state", "fingerprint_bitstring", "cluster_id", "cluster_label",
    "cluster_size", "cluster_medoid_flag", "distance_to_medoid",
    "dominant_key_interactions_json", "cluster_interpretation", "warnings_json",
]

FINAL_POSE_LABEL_COLUMNS = [
    "pose_id", "docking_task_id", "ligand_id", "target_receptor_id", "task_type",
    "docking_engine", "pose_rank", "rmsd_symmetry_corrected", "sanity_status",
    "stage3_pose_label", "ifp_tanimoto_to_native", "ifp_tanimoto_to_consensus",
    "key_interaction_recall_native", "key_interaction_recall_consensus",
    "interaction_recovery_label", "final_pose_label", "final_pose_label_confidence",
    "final_pose_label_reason", "warnings_json",
]

PLIP_CROSSCHECK_COLUMNS = [
    "analysis_id", "pose_id", "complex_id", "entity_type", "plip_status", "plip_version",
    "plip_interactions_json", "prolif_interactions_json", "agreement_score",
    "major_disagreements_json", "disagreement_reason", "use_for_report_flag",
    "warnings_json",
]

STAGE6_INTERACTION_FEATURE_COLUMNS = [
    "pose_id", "docking_task_id", "ligand_id", "target_receptor_id", "task_type",
    "docking_engine", "receptor_state", "pose_rank", "original_docking_score",
    "gnina_cnnscore", "gnina_cnnaffinity", "rmsd_symmetry_corrected", "sanity_status",
    "stage3_pose_label", "final_pose_label", "ifp_tanimoto_to_native",
    "ifp_tanimoto_to_consensus", "key_interaction_recall_native",
    "key_interaction_precision_native", "key_interaction_f1_native",
    "key_interaction_recall_consensus", "key_interaction_precision_consensus",
    "key_interaction_f1_consensus", "num_interactions", "num_key_interactions",
    "fingerprint_sparse_json", "binding_mode_cluster_id",
]
