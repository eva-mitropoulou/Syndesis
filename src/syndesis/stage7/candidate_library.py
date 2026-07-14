from __future__ import annotations

from pathlib import Path

import pandas as pd

from syndesis.common.io import write_json
from syndesis.stage7.library_export import (
    build_master_library,
    fetch_known_egfr_ligands,
    import_vendor_library,
    standardize_candidates,
)
from syndesis.stage7.load_stage_inputs import load_stage7_config, stage7_paths
from syndesis.stage7.report_stage7 import write_stage7_report


def clean_activity_data(config_path: str | Path) -> dict:
    return standardize_candidates(config_path)


def build_analog_series(config_path: str | Path) -> dict:
    config = load_stage7_config(config_path)
    paths = stage7_paths(config)
    aggregate = pd.read_parquet(paths["processed"] / "known_egfr_ligands_aggregated.parquet")
    if aggregate.empty:
        series = pd.DataFrame()
    else:
        series = aggregate.assign(
            series_id=aggregate["endpoint_type"].fillna("unknown") + "_" + aggregate["molecule_id"].str[:12],
            scaffold_id="",
            bemis_murcko_scaffold_smiles="",
            document_id="aggregate",
            assay_id="aggregate",
            activity_rank_within_series=1,
            series_size=1,
            source=aggregate["source_list"],
            include_in_retrospective_benchmark=False,
            warnings_json="[]",
        )[["series_id", "molecule_id", "scaffold_id", "bemis_murcko_scaffold_smiles", "endpoint_type", "assay_context_group", "document_id", "assay_id", "median_p_activity", "activity_rank_within_series", "series_size", "source", "include_in_retrospective_benchmark", "warnings_json"]]
    from syndesis.common.io import write_table

    write_table(paths["processed"] / "known_egfr_analog_series.parquet", series)
    write_table(paths["processed"] / "known_egfr_analog_series.csv", series)
    return {"status": "complete", "series_rows": int(len(series))}


def compute_candidate_similarity(config_path: str | Path) -> dict:
    return build_master_library(config_path)


def filter_candidate_library(config_path: str | Path) -> dict:
    summary = build_master_library(config_path)
    config = load_stage7_config(config_path)
    paths = stage7_paths(config)
    master = pd.read_parquet(paths["processed"] / "candidate_library_master.parquet")
    flags = master[["molecule_id", "standard_smiles", "source", "subsource", "mw", "clogp", "tpsa", "hbd", "hba", "rotatable_bonds", "formal_charge", "qed", "filter_pass", "hard_scope_pass", "soft_medchem_pass", "risk_flags_json", "include_in_screening_library"]].copy()
    for col in ["heavy_atom_count", "ring_count", "aromatic_ring_count", "fraction_sp3", "pains_flag", "brenk_flag", "reactive_flag", "aggregator_risk_flag", "covalent_warhead_flag", "egfr_cys797_warhead_flag", "allosteric_scope_flag", "macrocycle_flag", "metal_flag", "mixture_flag", "invalid_structure_flag", "property_window_pass", "exclusion_reason"]:
        if col not in flags:
            flags[col] = False if col.endswith("_flag") or col.endswith("_pass") else ""
    from syndesis.common.io import write_table

    write_table(paths["processed"] / "candidate_filter_flags.parquet", flags)
    write_table(paths["processed"] / "candidate_filter_flags.csv", flags)
    return summary


def select_screening_subsets(config_path: str | Path) -> dict:
    return build_master_library(config_path)


def prepare_candidate_ligands(config_path: str | Path) -> dict:
    return build_master_library(config_path)


def report_stage7(config_path: str | Path) -> dict:
    config = load_stage7_config(config_path)
    paths = stage7_paths(config)
    report = write_stage7_report(paths)
    payload = {"status": "complete", "report": str(report)}
    write_json(paths["processed"] / "stage7_report_summary.json", payload)
    return payload


def run_stage7_all(config_path: str | Path) -> dict:
    fetch_known_egfr_ligands(config_path)
    import_vendor_library(config_path)
    standardize_candidates(config_path)
    build_analog_series(config_path)
    filter_candidate_library(config_path)
    summary = report_stage7(config_path)
    return summary
