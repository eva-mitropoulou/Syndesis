from __future__ import annotations

from typing import Any

import pandas as pd


BENCHMARK_COLUMNS = [
    "complex_id",
    "pdb_id",
    "title",
    "primary_citation_doi",
    "release_date",
    "experimental_method",
    "resolution_angstrom",
    "r_work",
    "r_free",
    "organism",
    "uniprot_accession",
    "gene_name",
    "protein_name",
    "label_asym_id",
    "auth_asym_id",
    "entity_id",
    "chain_sequence_length",
    "modeled_residue_count",
    "deposited_residue_count",
    "kinase_domain_start_uniprot",
    "kinase_domain_end_uniprot",
    "kinase_domain_coverage_fraction",
    "mutation_flag",
    "mutation_list",
    "ligand_comp_id",
    "ligand_instance_id",
    "ligand_name",
    "ligand_formula",
    "ligand_inchi_key",
    "ligand_smiles",
    "ligand_heavy_atom_count",
    "ligand_mw",
    "ligand_occupancy_min",
    "ligand_occupancy_mean",
    "ligand_bfactor_mean",
    "pocket_bfactor_mean",
    "ligand_to_pocket_bfactor_ratio",
    "ligand_altloc_flag",
    "covalent_flag",
    "covalent_evidence",
    "warhead_flag",
    "warhead_reason",
    "atp_site_flag",
    "atp_site_evidence",
    "allosteric_flag",
    "ligand_class",
    "inhibitor_type_if_known",
    "klifs_structure_id",
    "klifs_kinase_id",
    "klifs_ligand_id",
    "klifs_pocket_id",
    "kincore_state",
    "dfg_state",
    "chelix_state",
    "saltbridge_state",
    "hrd_state",
    "activation_loop_state",
    "missing_active_site_residues",
    "missing_active_site_atoms",
    "active_site_completeness_score",
    "validation_report_available",
    "ligand_validation_available",
    "validation_xml_path",
    "validation_pdf_path",
    "raw_mmcif_path",
    "native_complex_path",
    "native_ligand_sdf_path",
    "receptor_clean_path",
    "pocket_waters_path",
    "quality_tier",
    "quality_score",
    "include_in_stage1_benchmark",
    "exclusion_reason",
    "warnings_json",
]

QUALITY_TIERS = {"Tier A", "Tier B", "Tier C", "Rejected"}


def empty_benchmark_frame() -> pd.DataFrame:
    return pd.DataFrame(columns=BENCHMARK_COLUMNS)


def normalize_record(record: dict[str, Any]) -> dict[str, Any]:
    return {column: record.get(column) for column in BENCHMARK_COLUMNS}


def benchmark_frame(records: list[dict[str, Any]]) -> pd.DataFrame:
    return pd.DataFrame([normalize_record(record) for record in records], columns=BENCHMARK_COLUMNS)

