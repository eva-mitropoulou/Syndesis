from __future__ import annotations

import json
from datetime import date
from pathlib import Path
from typing import Any

import pandas as pd
from rdkit import Chem

from egfr_dockingforge.common.io import write_json, write_table
from egfr_dockingforge.stage7.activity_aggregation import aggregate_known
from egfr_dockingforge.stage7.activity_cleaning import p_activity_from_nm, to_nm
from egfr_dockingforge.stage7.chembl_client import discover_egfr_target, fetch_activities
from egfr_dockingforge.stage7.diversity_selection import select_subsets
from egfr_dockingforge.stage7.filtering import flags
from egfr_dockingforge.stage7.ligand_preparation import prepare_ligands
from egfr_dockingforge.stage7.load_stage_inputs import load_stage7_config, load_stage7_inputs, stage7_paths
from egfr_dockingforge.stage7.scaffold_analysis import murcko
from egfr_dockingforge.stage7.similarity import closest, novelty_bucket, prepare_known
from egfr_dockingforge.stage7.standardize import molecule_id_from_inchikey, standardize_smiles
from egfr_dockingforge.stage7.vendor_import import fetch_cartblanche_similarity, read_vendor_file


def _load(config_path: str | Path) -> tuple[dict[str, Any], dict[str, Path], dict[str, pd.DataFrame]]:
    config = load_stage7_config(config_path)
    paths = stage7_paths(config)
    inputs = load_stage7_inputs(config)
    return config, paths, inputs


def fetch_known_egfr_ligands(config_path: str | Path) -> dict[str, Any]:
    config, paths, inputs = _load(config_path)
    target = discover_egfr_target(config["known_ligands"]["uniprot_id"])
    chembl = pd.DataFrame(fetch_activities(target, int(config["known_ligands"]["chembl_activity_limit"])))
    chembl["retrieval_date"] = date.today().isoformat()
    write_table(paths["raw"] / "chembl" / "chembl_egfr_activities_raw.parquet", chembl)
    write_table(paths["raw"] / "chembl" / "chembl_egfr_activities_raw.csv", chembl)
    native = _stage1_native_ligands(inputs["cocrystal_benchmark"])
    write_table(paths["raw"] / "stage1_native_ligands_raw.parquet", native)
    bindingdb = pd.DataFrame()
    try:
        from egfr_dockingforge.stage7.bindingdb_client import fetch_bindingdb_uniprot

        bindingdb = pd.DataFrame(fetch_bindingdb_uniprot(config["known_ligands"]["uniprot_id"], int(config["known_ligands"]["bindingdb_cutoff_nm"]), int(config["known_ligands"]["bindingdb_timeout_seconds"]), int(config["known_ligands"]["bindingdb_max_bytes"])))
    except Exception as exc:
        bindingdb = pd.DataFrame([{"fetch_status": "failed", "error": str(exc)}])
    write_table(paths["raw"] / "bindingdb" / "bindingdb_egfr_activities_raw.parquet", bindingdb)
    write_table(paths["raw"] / "bindingdb" / "bindingdb_egfr_activities_raw.csv", bindingdb)
    return {"status": "complete", "chembl_rows": int(len(chembl)), "bindingdb_rows": int(len(bindingdb)), "native_rows": int(len(native))}


def import_vendor_library(config_path: str | Path) -> dict[str, Any]:
    config, paths, _inputs = _load(config_path)
    manifest_rows = []
    rows = []
    for file_path in config["vendor"].get("local_files", []):
        frame = read_vendor_file(file_path)
        rows.extend(_vendor_rows(frame, "local_vendor_file", str(file_path)))
        manifest_rows.append(_manifest("other_vendor", str(file_path), len(frame), "success", []))
    if config["vendor"].get("cartblanche_enabled", False):
        for smiles in config["vendor"].get("cartblanche_query_smiles", []):
            try:
                frame = fetch_cartblanche_similarity(smiles, int(config["vendor"]["cartblanche_distance"]))
                rows.extend(_vendor_rows(frame, "cartblanche_similarity", smiles))
                manifest_rows.append(_manifest("zinc_vendor", smiles, len(frame), "success", []))
            except Exception as exc:
                manifest_rows.append(_manifest("zinc_vendor", smiles, 0, "failed", [str(exc)]))
    vendor = pd.DataFrame(rows).head(int(config["vendor"]["max_vendor_records"]))
    manifest = pd.DataFrame(manifest_rows)
    write_table(paths["raw"] / "vendor" / "vendor_manifest.parquet", manifest)
    write_table(paths["raw"] / "vendor" / "vendor_manifest.csv", manifest)
    write_table(paths["processed"] / "vendor_molecules_standardized.parquet", vendor)
    write_table(paths["processed"] / "vendor_molecules_standardized.csv", vendor)
    return {"status": "complete", "vendor_rows": int(len(vendor)), "manifest_rows": int(len(manifest))}


def standardize_candidates(config_path: str | Path) -> dict[str, Any]:
    config, paths, _inputs = _load(config_path)
    chembl_raw = pd.read_parquet(paths["raw"] / "chembl" / "chembl_egfr_activities_raw.parquet")
    native_raw = pd.read_parquet(paths["raw"] / "stage1_native_ligands_raw.parquet")
    std_rows = []
    meas_rows = []
    for row in chembl_raw.to_dict("records"):
        raw_id = str(row.get("activity_id") or row.get("activity_chembl_id") or row.get("record_id"))
        smi = row.get("canonical_smiles") or row.get("molecule_structures", {}).get("canonical_smiles") if isinstance(row.get("molecule_structures"), dict) else row.get("canonical_smiles")
        std = standardize_smiles(smi, raw_id, "chembl_known_ligand")
        std_rows.append(std)
        if std["standardization_status"] == "success":
            value_nm = to_nm(row.get("standard_value"), row.get("standard_units"))
            meas_rows.append(_measurement(row, std, value_nm, "chembl_known_ligand"))
    for row in native_raw.to_dict("records"):
        std = standardize_smiles(row["ligand_smiles"], row["complex_id"], "stage1_native_ligand")
        std_rows.append(std)
    stdlog = pd.DataFrame(std_rows)
    measurements = pd.DataFrame(meas_rows)
    write_table(paths["processed"] / "molecule_standardization_log.parquet", stdlog)
    write_table(paths["processed"] / "molecule_standardization_log.csv", stdlog)
    write_table(paths["processed"] / "known_egfr_activity_measurements.parquet", measurements)
    write_table(paths["processed"] / "known_egfr_activity_measurements.csv", measurements)
    aggregate = aggregate_known(measurements)
    write_table(paths["processed"] / "known_egfr_ligands_aggregated.parquet", aggregate)
    write_table(paths["processed"] / "known_egfr_ligands_aggregated.csv", aggregate)
    return {"status": "complete", "standardized_rows": int(len(stdlog)), "measurement_rows": int(len(measurements))}


def build_master_library(config_path: str | Path) -> dict[str, Any]:
    config, paths, inputs = _load(config_path)
    if not (paths["processed"] / "known_egfr_ligands_aggregated.parquet").exists():
        standardize_candidates(config_path)
    aggregate = pd.read_parquet(paths["processed"] / "known_egfr_ligands_aggregated.parquet")
    native = _stage1_native_ligands(inputs["cocrystal_benchmark"])
    vendor = pd.read_parquet(paths["processed"] / "vendor_molecules_standardized.parquet") if (paths["processed"] / "vendor_molecules_standardized.parquet").exists() else pd.DataFrame()
    rows = []
    for row in aggregate.to_dict("records"):
        rows.append(_master_row(row["molecule_id"], "chembl_known_ligand", "chembl", "known_activity_reference", row["standard_smiles"], row["inchi_key"], config, row))
    for row in native.to_dict("records"):
        std = standardize_smiles(row["ligand_smiles"], row["complex_id"], "stage1_native_ligand")
        if std["standardization_status"] == "success":
            rows.append(_master_row(molecule_id_from_inchikey(std["inchi_key"], std["standard_smiles"]), "stage1_native_ligand", row["pdb_id"], "native_pose_reference", std["standard_smiles"], std["inchi_key"], config, {"native_ligand_flag": True}))
    for row in vendor.to_dict("records"):
        if row.get("standardization_status") == "success":
            rows.append(_master_row(row["molecule_id"], row["source"], row["subsource"], "prospective_candidate", row["standard_smiles"], row["inchi_key"], config, {"vendor_availability_status": row.get("availability_status", "unknown")}))
    master = pd.DataFrame(rows).drop_duplicates("molecule_id")
    known_refs = [
        (r["molecule_id"], r["source"], r["standard_smiles"], r.get("median_standard_value_nM"), r.get("endpoint_type_if_known"))
        for r in master[master["source"].isin(["chembl_known_ligand", "stage1_native_ligand"])].to_dict("records")
    ]
    known_refs_prepared = prepare_known(known_refs)
    known_scaffolds = set(master[master["source"].isin(["chembl_known_ligand", "stage1_native_ligand"])]["bemis_murcko_scaffold_smiles"])
    similarity_rows = []
    for idx, row in master.iterrows():
        c_id, c_source, c_act, c_endpoint, sim = closest(row["standard_smiles"], known_refs_prepared)
        same_key = bool(row["inchi_key"] and row["inchi_key"] in set(master[master["source"].isin(["chembl_known_ligand", "stage1_native_ligand"])]["inchi_key"]))
        scaffold_match = row["bemis_murcko_scaffold_smiles"] in known_scaffolds
        bucket = novelty_bucket(sim, same_key and row["source"] not in {"chembl_known_ligand", "stage1_native_ligand"}, scaffold_match)
        master.loc[idx, "closest_known_molecule_id"] = c_id
        master.loc[idx, "tanimoto_to_closest_known"] = sim
        master.loc[idx, "novelty_bucket"] = bucket
        similarity_rows.append({"molecule_id": row["molecule_id"], "source": row["source"], "standard_smiles": row["standard_smiles"], "scaffold_id": row["scaffold_id"], "bemis_murcko_scaffold_smiles": row["bemis_murcko_scaffold_smiles"], "closest_known_molecule_id": c_id, "closest_known_source": c_source, "closest_known_activity_nM": c_act, "closest_known_endpoint_type": c_endpoint, "tanimoto_to_closest_known": sim, "scaffold_match_to_known_flag": scaffold_match, "native_ligand_similarity_max": sim, "known_active_similarity_max": sim, "novelty_bucket": bucket, "hard_novel_flag": bucket == "scaffold_novel", "analog_like_flag": bucket == "close_analog", "potential_contamination_flag": bucket == "known_duplicate", "warnings_json": json.dumps([])})
    sim_df = pd.DataFrame(similarity_rows)
    write_table(paths["processed"] / "candidate_similarity_to_known.parquet", sim_df)
    write_table(paths["processed"] / "candidate_similarity_to_known.csv", sim_df)
    write_table(paths["processed"] / "candidate_library_master.parquet", master)
    write_table(paths["processed"] / "candidate_library_master.csv", master)
    subsets = select_subsets(master, config)
    write_table(paths["processed"] / "screening_subsets.parquet", subsets)
    write_table(paths["processed"] / "screening_subsets.csv", subsets)
    prep = prepare_ligands(master, subsets, config, paths)
    write_table(paths["processed"] / "candidate_ligand_preparation.parquet", prep)
    write_table(paths["processed"] / "candidate_ligand_preparation.csv", prep)
    stage8 = _stage8_input(prep, master)
    write_table(paths["processed"] / "stage8_screening_input.parquet", stage8)
    write_table(paths["processed"] / "stage8_screening_input.csv", stage8)
    write_json(paths["processed"] / "stage7_summary.json", {"status": "complete", "master_rows": int(len(master)), "stage8_rows": int(len(stage8))})
    return {"status": "complete", "master_rows": int(len(master)), "stage8_rows": int(len(stage8))}


def _stage1_native_ligands(benchmark: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for row in benchmark[benchmark["include_in_stage1_benchmark"].astype(bool)].to_dict("records"):
        smiles = row.get("ligand_smiles") or _smiles_from_file(row.get("native_ligand_sdf_path"))
        rows.append({"complex_id": row["complex_id"], "pdb_id": row["pdb_id"], "source": "stage1_native_ligand", "ligand_smiles": smiles, "inchi_key": row.get("ligand_inchi_key") or ""})
    return pd.DataFrame(rows)


def _smiles_from_file(path: str | None) -> str:
    if not path:
        return ""
    mol = Chem.MolFromPDBFile(path, sanitize=True, removeHs=True)
    if mol is None:
        return ""
    return Chem.MolToSmiles(mol, canonical=True)


def _vendor_rows(frame: pd.DataFrame, subsource: str, origin: str) -> list[dict]:
    rows = []
    for idx, row in frame.iterrows():
        smi = row.get("smiles") or row.get("SMILES")
        std = standardize_smiles(smi, f"vendor_{idx}", "zinc_vendor")
        if std["standardization_status"] == "success":
            rows.append({"vendor_molecule_id": row.get("zinc_id", f"vendor_{idx}"), "source": "zinc_vendor", "subsource": subsource, "vendor_catalog_id": row.get("zinc_id", ""), "availability_status": "cartblanche_reported", "purchasable_flag": True, "make_on_demand_flag": True, "supplier_name": str(row.get("catalogs", "")), "smiles_raw": smi, "standard_smiles": std["standard_smiles"], "inchi_key": std["inchi_key"], "molecule_id": molecule_id_from_inchikey(std["inchi_key"], std["standard_smiles"]), "standardization_status": "success", "import_warnings_json": std["warnings_json"], "source_file_or_query": origin})
    return rows


def _manifest(source: str, origin: str, rows: int, status: str, warnings: list[str]) -> dict:
    return {"vendor_source": source, "vendor_file": origin, "source_version": "ZINC22/CartBlanche22", "download_url_or_origin": origin, "retrieval_date": date.today().isoformat(), "license_or_terms_note": "ZINC academic/research terms; verify procurement terms before purchase.", "num_records_raw": rows, "checksum": "", "import_status": status, "warnings_json": json.dumps(warnings)}


def _measurement(row: dict, std: dict, value_nm: float | None, source: str) -> dict:
    mol_id = molecule_id_from_inchikey(std["inchi_key"], std["standard_smiles"])
    endpoint = row.get("standard_type")
    return {"measurement_id": str(row.get("activity_id")), "molecule_id": mol_id, "source": source, "source_record_id": str(row.get("activity_id")), "source_version": "ChEMBL API current", "target_name": "EGFR", "target_uniprot": "P00533", "target_chembl_id": row.get("target_chembl_id"), "bindingdb_target_id": "", "organism": row.get("target_organism"), "mutation_status": "", "mutation_list": "", "assay_id": row.get("assay_chembl_id"), "assay_description": row.get("assay_description"), "assay_type": row.get("assay_type"), "assay_cell_type": "", "assay_organism": "", "endpoint_type": endpoint, "relation": row.get("standard_relation"), "value_original": row.get("standard_value"), "unit_original": row.get("standard_units"), "standard_value_nM": value_nm, "p_activity": p_activity_from_nm(value_nm), "activity_direction": "lower_is_stronger", "document_id": row.get("document_chembl_id"), "doi": row.get("document_journal"), "pubmed_id": "", "patent_id": "", "molecule_smiles_raw": std["raw_smiles"], "molecule_smiles_standard": std["standard_smiles"], "inchi_key": std["inchi_key"], "standardization_status": std["standardization_status"], "include_in_known_activity_set": value_nm is not None and row.get("standard_relation") == "=", "exclusion_reason": "" if value_nm is not None else "missing_standard_value", "warnings_json": json.dumps([])}


def _master_row(molecule_id: str, source: str, subsource: str, role: str, smiles: str, inchi_key: str, config: dict, extra: dict) -> dict:
    filt = flags(smiles, source, config)
    scaffold = murcko(smiles)
    return {"molecule_id": molecule_id, "source": source, "subsource": subsource, "screening_role": role, "source_record_ids_json": extra.get("source_record_ids_json", "[]"), "source_version": "stage7", "standard_smiles": smiles, "inchi_key": inchi_key, "inchi_key_connectivity": inchi_key.split("-")[0] if inchi_key else "", "scaffold_id": scaffold, "bemis_murcko_scaffold_smiles": scaffold, "known_activity_status": extra.get("known_activity_status", "unknown"), "median_p_activity_if_known": extra.get("median_p_activity"), "endpoint_type_if_known": extra.get("endpoint_type"), "native_ligand_flag": bool(extra.get("native_ligand_flag", False)), "approved_or_clinical_flag": bool(extra.get("approved_or_clinical_flag", False)), "vendor_availability_status": extra.get("vendor_availability_status", ""), "parent_molecule_id": "", "generation_method": "", "filter_pass": filt["include_in_screening_library"], "soft_medchem_pass": filt["soft_medchem_pass"], "closest_known_molecule_id": "", "tanimoto_to_closest_known": None, "novelty_bucket": "", "screening_subset_list_json": "[]", "prepared_ligand_available_flag": False, "warnings_json": "[]", **{k: filt.get(k) for k in ["mw", "clogp", "tpsa", "hbd", "hba", "rotatable_bonds", "formal_charge", "qed", "hard_scope_pass", "risk_flags_json", "include_in_screening_library"]}}


def _stage8_input(prep: pd.DataFrame, master: pd.DataFrame) -> pd.DataFrame:
    if prep.empty:
        return pd.DataFrame()
    merged = prep.merge(master[["molecule_id", "screening_role", "mw", "clogp", "tpsa", "novelty_bucket", "closest_known_molecule_id", "tanimoto_to_closest_known", "risk_flags_json"]], on="molecule_id", how="left")
    return merged.rename(columns={"pdbqt_path_if_available": "docking_format_path_if_available"})[["prepared_ligand_id", "molecule_id", "source", "screening_role", "screening_subset", "standard_smiles", "prepared_smiles", "sdf_path", "docking_format_path_if_available", "protonation_state_id", "tautomer_state_id", "conformer_id", "mw", "clogp", "tpsa", "novelty_bucket", "closest_known_molecule_id", "tanimoto_to_closest_known", "risk_flags_json"]]
