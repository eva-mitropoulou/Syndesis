from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pandas as pd
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from joblib import Parallel, delayed

from egfr_dockingforge.common.io import ensure_dir, load_yaml, project_root, resolve_path, write_json, write_table
from egfr_dockingforge.stage1.complex_classification import classify_complex
from egfr_dockingforge.stage1.file_export import export_reference_files
from egfr_dockingforge.stage1.kincore_metadata import infer_kincore_metadata
from egfr_dockingforge.stage1.klifs_client import fetch_klifs_metadata
from egfr_dockingforge.stage1.ligand_extraction import find_ligands_for_chain
from egfr_dockingforge.stage1.mmcif_parser import (
    chem_comp_descriptors,
    chem_comp_table,
    entry_methods,
    entry_resolution,
    parse_mmcif,
    protein_chains,
    residue_count,
    r_factor,
)
from egfr_dockingforge.stage1.quality_scoring import assign_quality
from egfr_dockingforge.stage1.rcsb_client import (
    fetch_chem_comp,
    fetch_entry_metadata,
    fetch_mmcif,
    fetch_validation_files,
    file_record,
    search_pdb_ids,
    write_download_manifest,
)
from egfr_dockingforge.stage1.report_stage1 import write_stage1_report
from egfr_dockingforge.stage1.schemas import benchmark_frame, empty_benchmark_frame


def config_paths(config: dict[str, Any], root: Path) -> dict[str, Path]:
    return {key: resolve_path(value, root) for key, value in config["paths"].items()}


def load_stage1_config(config_path: str | Path) -> dict[str, Any]:
    return load_yaml(resolve_path(config_path, project_root()))


def metadata_value(metadata: dict[str, Any], path: list[str], default: Any = None) -> Any:
    value: Any = metadata
    for key in path:
        if not isinstance(value, dict) or key not in value:
            return default
        value = value[key]
    return value


def first_list_value(value: Any) -> Any:
    if isinstance(value, list):
        return value[0] if value else None
    return value


def entry_title(metadata: dict[str, Any]) -> str | None:
    return metadata_value(metadata, ["struct", "title"])


def entry_doi(metadata: dict[str, Any]) -> str | None:
    citations = metadata.get("citation") or []
    if isinstance(citations, list):
        for citation in citations:
            doi = citation.get("pdbx_database_id_DOI") if isinstance(citation, dict) else None
            if doi not in {None, "", "?", "."}:
                return str(doi)
    return None


def entry_release_date(metadata: dict[str, Any]) -> str | None:
    return metadata_value(metadata, ["rcsb_accession_info", "initial_release_date"])


def metadata_resolution(metadata: dict[str, Any]) -> float | None:
    values = metadata_value(metadata, ["rcsb_entry_info", "resolution_combined"], [])
    value = first_list_value(values)
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def metadata_method(metadata: dict[str, Any]) -> str | None:
    methods = metadata_value(metadata, ["exptl"], [])
    if isinstance(methods, list):
        names = [row.get("method") for row in methods if isinstance(row, dict) and row.get("method")]
        return "; ".join(names) if names else None
    return None


def chain_metadata(chain: Any, config: dict[str, Any], pdb_id: str) -> dict[str, Any]:
    count = residue_count(chain)
    target = config["target"]
    mutations = config.get("known_mutations", {}).get(pdb_id.upper(), [])
    return {
        "organism": target.get("organism"),
        "uniprot_accession": target.get("uniprot_accession"),
        "gene_name": "EGFR",
        "protein_name": "Epidermal growth factor receptor",
        "entity_id": None,
        "chain_sequence_length": count,
        "modeled_residue_count": count,
        "deposited_residue_count": count,
        "kinase_domain_start_uniprot": target["kinase_domain_uniprot"]["start"],
        "kinase_domain_end_uniprot": target["kinase_domain_uniprot"]["end"],
        "kinase_domain_coverage_fraction": min(1.0, count / 268.0) if count else 0.0,
        "mutation_flag": bool(mutations),
        "mutation_list": ",".join(mutations) if mutations else "",
    }


def complex_id_for(pdb_id: str, chain_id: str, ligand_comp_id: str, ligand_instance_id: str) -> str:
    return f"{pdb_id.upper()}_{chain_id}_{ligand_comp_id}_{ligand_instance_id}"


def ligand_ccd_metadata(comp_id: str, config: dict[str, Any], ccd_dir: Path) -> dict[str, Any]:
    ccd_path = fetch_chem_comp(comp_id, config, ccd_dir)
    if not ccd_path:
        return {}
    mmcif = MMCIF2Dict(str(ccd_path))
    descriptors = chem_comp_descriptors(mmcif).get(comp_id.upper(), {})
    table = chem_comp_table(mmcif).get(comp_id.upper(), {})
    return {**table, **descriptors}


def acquire_stage1(config_path: str | Path) -> dict[str, Any]:
    root = project_root()
    config = load_stage1_config(config_path)
    paths = config_paths(config, root)
    for key in ("raw_rcsb_mmcif", "raw_rcsb_metadata", "raw_rcsb_validation", "raw_klifs", "raw_ligands"):
        ensure_dir(paths[key])
    pdb_ids = search_pdb_ids(config)
    force = bool(config["rcsb"].get("downloads", {}).get("force", False))
    records: list[dict[str, Any]] = []
    for pdb_id in pdb_ids:
        mmcif_path = fetch_mmcif(pdb_id, config, paths["raw_rcsb_mmcif"], force=force)
        metadata = fetch_entry_metadata(pdb_id, config, paths["raw_rcsb_metadata"], force=force)
        xml_path, pdf_path = fetch_validation_files(pdb_id, config, paths["raw_rcsb_validation"], force=force)
        fetch_klifs_metadata(pdb_id, config, paths["raw_klifs"])
        parser_version = config["stage"]["parser_version"]
        records.append(
            {
                "pdb_id": pdb_id,
                "mmcif": file_record(mmcif_path, config["rcsb"]["mmcif_url_template"].format(pdb_id=pdb_id), parser_version),
                "metadata_path": str(paths["raw_rcsb_metadata"] / f"{pdb_id.lower()}_entry.json"),
                "validation_xml_path": str(xml_path) if xml_path else None,
                "validation_pdf_path": str(pdf_path) if pdf_path else None,
                "metadata_title": entry_title(metadata),
            }
        )
    manifest_path = paths["interim"] / "stage1_acquisition_manifest.json"
    write_download_manifest(records, manifest_path)
    return {"pdb_ids": pdb_ids, "manifest_path": str(manifest_path)}


def process_pdb(pdb_id: str, config: dict[str, Any], paths: dict[str, Path]) -> list[dict[str, Any]]:
    warnings: list[str] = []
    mmcif_path = paths["raw_rcsb_mmcif"] / f"{pdb_id.lower()}.cif"
    metadata_path = paths["raw_rcsb_metadata"] / f"{pdb_id.lower()}_entry.json"
    if not mmcif_path.exists():
        fetch_mmcif(pdb_id, config, paths["raw_rcsb_mmcif"], force=False)
    if not metadata_path.exists():
        metadata = fetch_entry_metadata(pdb_id, config, paths["raw_rcsb_metadata"], force=False)
    else:
        metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
    xml_path = paths["raw_rcsb_validation"] / f"{pdb_id.lower()}_validation.xml"
    pdf_path = paths["raw_rcsb_validation"] / f"{pdb_id.lower()}_validation.pdf"
    parsed = parse_mmcif(mmcif_path)
    descriptors = chem_comp_descriptors(parsed.mmcif)
    compounds = chem_comp_table(parsed.mmcif)
    klifs, klifs_warnings = fetch_klifs_metadata(pdb_id, config, paths["raw_klifs"])
    warnings.extend(klifs_warnings)

    resolution = metadata_resolution(metadata) or entry_resolution(parsed.mmcif)
    method = metadata_method(metadata) or entry_methods(parsed.mmcif)
    title = entry_title(metadata)
    records: list[dict[str, Any]] = []
    for chain in protein_chains(parsed.structure):
        chain_meta = chain_metadata(chain, config, pdb_id)
        ligands = find_ligands_for_chain(chain, config)
        for ligand in ligands:
            ligand_meta = compounds.get(ligand.comp_id, {})
            ligand_desc = descriptors.get(ligand.comp_id, {})
            ccd_meta = ligand_ccd_metadata(ligand.comp_id, config, paths["raw_ligands"])
            ligand_meta = {**ligand_meta, **ccd_meta}
            ligand_desc = {**ligand_desc, **{k: ccd_meta.get(k) for k in ("smiles", "inchi_key") if ccd_meta.get(k)}}
            cid = complex_id_for(pdb_id, chain.id, ligand.comp_id, ligand.instance_id)
            kincore, kincore_warnings = infer_kincore_metadata(pdb_id, chain.id, config, klifs)
            row_warnings = list(warnings) + kincore_warnings
            bfactor_ratio = None
            if ligand.bfactor_mean is not None and ligand.pocket_bfactor_mean:
                bfactor_ratio = ligand.bfactor_mean / ligand.pocket_bfactor_mean
            classification = classify_complex(
                pdb_id=pdb_id,
                ligand=ligand,
                chain=chain,
                mmcif=parsed.mmcif,
                ligand_smiles=ligand_desc.get("smiles"),
                config=config,
            )
            exports = export_reference_files(
                parsed.structure,
                cid,
                chain.id,
                ligand.residue_id,
                paths["reference_complexes"],
            )
            record: dict[str, Any] = {
                "complex_id": cid,
                "pdb_id": pdb_id,
                "title": title,
                "primary_citation_doi": entry_doi(metadata),
                "release_date": entry_release_date(metadata),
                "experimental_method": method,
                "resolution_angstrom": resolution,
                "r_work": r_factor(parsed.mmcif, "_refine.ls_R_factor_R_work"),
                "r_free": r_factor(parsed.mmcif, "_refine.ls_R_factor_R_free"),
                "label_asym_id": ligand.label_asym_id,
                "auth_asym_id": chain.id,
                "ligand_comp_id": ligand.comp_id,
                "ligand_instance_id": ligand.instance_id,
                "ligand_name": ligand_meta.get("ligand_name"),
                "ligand_formula": ligand_meta.get("ligand_formula"),
                "ligand_inchi_key": ligand_desc.get("inchi_key"),
                "ligand_smiles": ligand_desc.get("smiles"),
                "ligand_heavy_atom_count": ligand.heavy_atom_count,
                "ligand_mw": ligand_meta.get("ligand_mw"),
                "ligand_occupancy_min": ligand.occupancy_min,
                "ligand_occupancy_mean": ligand.occupancy_mean,
                "ligand_bfactor_mean": ligand.bfactor_mean,
                "pocket_bfactor_mean": ligand.pocket_bfactor_mean,
                "ligand_to_pocket_bfactor_ratio": bfactor_ratio,
                "ligand_altloc_flag": ligand.altloc_flag,
                "validation_report_available": xml_path.exists() or pdf_path.exists(),
                "ligand_validation_available": xml_path.exists(),
                "validation_xml_path": str(xml_path) if xml_path.exists() else None,
                "validation_pdf_path": str(pdf_path) if pdf_path.exists() else None,
                "raw_mmcif_path": str(mmcif_path),
                "warnings_json": None,
                **chain_meta,
                **classification,
                **klifs,
                **kincore,
                **exports,
            }
            quality = assign_quality(record, config)
            record.update(quality)
            row_warnings.extend(classification.get("classification_warnings") or [])
            if not ligand_desc.get("smiles"):
                row_warnings.append(f"No SMILES available for ligand {ligand.comp_id}.")
            if not record.get("klifs_structure_id"):
                row_warnings.append("KLIFS metadata unavailable; nullable KLIFS fields retained.")
            record["warnings_json"] = json.dumps(row_warnings, sort_keys=True)
            record.pop("hard_exclusion_reasons", None)
            record.pop("classification_warnings", None)
            records.append(record)
    if not records:
        records.append(
            {
                "complex_id": f"{pdb_id}_NO_USABLE_COMPLEX",
                "pdb_id": pdb_id,
                "quality_tier": "Rejected",
                "include_in_stage1_benchmark": False,
                "exclusion_reason": "No usable EGFR chain/native ligand complex after parsing.",
                "warnings_json": json.dumps(warnings),
                "raw_mmcif_path": str(mmcif_path),
                "title": title,
                "experimental_method": method,
                "resolution_angstrom": resolution,
            }
        )
    return records


def write_stage1_tables(frame: pd.DataFrame, paths: dict[str, Path]) -> dict[str, str]:
    ensure_dir(paths["processed"])
    ensure_dir(paths["interim"])
    bool_columns = [
        "mutation_flag",
        "ligand_altloc_flag",
        "covalent_flag",
        "warhead_flag",
        "atp_site_flag",
        "allosteric_flag",
        "validation_report_available",
        "ligand_validation_available",
        "include_in_stage1_benchmark",
    ]
    for column in bool_columns:
        if column in frame.columns:
            frame[column] = frame[column].astype("boolean")
    benchmark_path = paths["processed"] / "egfr_cocrystal_benchmark.parquet"
    benchmark_csv = paths["processed"] / "egfr_cocrystal_benchmark.csv"
    rejected = frame[frame["include_in_stage1_benchmark"] != True].copy()
    rejected_path = paths["interim"] / "rejected_complexes.parquet"
    rejected_csv = paths["interim"] / "rejected_complexes.csv"
    write_table(benchmark_csv, frame)
    write_table(rejected_csv, rejected)
    try:
        write_table(benchmark_path, frame)
        write_table(rejected_path, rejected)
    except RuntimeError:
        # CSV sidecars are intentionally written first for environments without parquet engines.
        pass
    return {
        "benchmark_parquet": str(benchmark_path),
        "benchmark_csv": str(benchmark_csv),
        "rejected_parquet": str(rejected_path),
        "rejected_csv": str(rejected_csv),
    }


def build_stage1_benchmark(config_path: str | Path) -> dict[str, Any]:
    root = project_root()
    config = load_stage1_config(config_path)
    paths = config_paths(config, root)
    for path in paths.values():
        ensure_dir(path)
    acquire = acquire_stage1(config_path)
    pdb_ids = acquire["pdb_ids"]
    workers = int(config.get("runtime", {}).get("parse_workers", 8))
    nested = Parallel(n_jobs=max(1, workers), prefer="threads")(
        delayed(process_pdb)(pdb_id, config, paths) for pdb_id in pdb_ids
    )
    records = [record for rows in nested for record in rows]
    frame = benchmark_frame(records) if records else empty_benchmark_frame()
    outputs = write_stage1_tables(frame, paths)
    rejected = frame[frame["include_in_stage1_benchmark"] != True].copy() if not frame.empty else frame
    report_path = write_stage1_report(frame, rejected, paths["reports"] / "01_cocrystal_benchmark.html")
    summary = {
        "status": "complete",
        "candidate_pdb_entries": len(pdb_ids),
        "candidate_complexes": int(len(frame)),
        "retained_complexes": int((frame["include_in_stage1_benchmark"] == True).sum()) if not frame.empty else 0,
        "rejected_complexes": int(len(rejected)),
        "outputs": outputs,
        "report": str(report_path),
    }
    write_json(paths["processed"] / "stage1_summary.json", summary)
    return summary


def extract_native_ligands(config_path: str | Path) -> dict[str, Any]:
    return build_stage1_benchmark(config_path)


def export_reference_complexes(config_path: str | Path) -> dict[str, Any]:
    return build_stage1_benchmark(config_path)


def report_stage1(config_path: str | Path) -> dict[str, Any]:
    root = project_root()
    config = load_stage1_config(config_path)
    paths = config_paths(config, root)
    csv_path = paths["processed"] / "egfr_cocrystal_benchmark.csv"
    if not csv_path.exists():
        return build_stage1_benchmark(config_path)
    frame = pd.read_csv(csv_path)
    rejected = frame[frame["include_in_stage1_benchmark"] != True].copy()
    report = write_stage1_report(frame, rejected, paths["reports"] / "01_cocrystal_benchmark.html")
    return {"status": "complete", "report": str(report)}
