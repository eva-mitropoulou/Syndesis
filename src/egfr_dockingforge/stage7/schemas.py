from __future__ import annotations

MASTER_COLUMNS = [
    "molecule_id", "source", "subsource", "screening_role", "source_record_ids_json",
    "source_version", "standard_smiles", "inchi_key", "inchi_key_connectivity",
    "scaffold_id", "bemis_murcko_scaffold_smiles", "known_activity_status",
    "median_p_activity_if_known", "endpoint_type_if_known", "native_ligand_flag",
    "approved_or_clinical_flag", "vendor_availability_status", "parent_molecule_id",
    "generation_method", "mw", "clogp", "tpsa", "hbd", "hba", "rotatable_bonds",
    "formal_charge", "qed", "filter_pass", "hard_scope_pass", "soft_medchem_pass",
    "risk_flags_json", "closest_known_molecule_id", "tanimoto_to_closest_known",
    "novelty_bucket", "include_in_screening_library", "screening_subset_list_json",
    "prepared_ligand_available_flag", "warnings_json",
]

STAGE8_COLUMNS = [
    "prepared_ligand_id", "molecule_id", "source", "screening_role", "screening_subset",
    "standard_smiles", "prepared_smiles", "sdf_path", "docking_format_path_if_available",
    "protonation_state_id", "tautomer_state_id", "conformer_id", "mw", "clogp", "tpsa",
    "novelty_bucket", "closest_known_molecule_id", "tanimoto_to_closest_known",
    "risk_flags_json",
]

ALLOWED_SOURCES = {
    "stage1_native_ligand", "chembl_known_ligand", "bindingdb_known_ligand",
    "pubchem_known_ligand", "iuphar_or_clinical_control", "zinc_vendor",
    "enamine_vendor", "other_vendor", "generated_analog", "manual_analog",
}
