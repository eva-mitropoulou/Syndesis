from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pandas as pd

from syndesis.common.io import write_table
from syndesis.stage12.load_stage_inputs import load_stage12_config, load_stage12_inputs, stage12_paths
from syndesis.stage12.nonclaim_generator import nonclaim_text
from syndesis.stage12.risk_summary import merge_warnings
from syndesis.stage12.schemas import FINAL_CANDIDATE_COLUMNS, FINAL_RANKED_COLUMNS, validate_final_candidates


def _value(row: pd.Series, name: str, default: Any = None) -> Any:
    return row[name] if name in row and pd.notna(row[name]) else default


def _md_lookup(inputs: dict[str, pd.DataFrame | None]) -> dict[str, dict[str, Any]]:
    labels = inputs.get("stage11_md_stability")
    metrics = inputs.get("stage11_md_metrics")
    parameterization = inputs.get("stage11_parameterization_report")
    out: dict[str, dict[str, Any]] = {}
    if labels is not None:
        for row in labels.to_dict("records"):
            out.setdefault(row["molecule_id"], {}).update(row)
    if (
        metrics is not None
        and labels is not None
        and {"molecule_id", "md_candidate_id"} <= set(labels.columns)
        and "md_candidate_id" in metrics.columns
    ):
        # Guard the columns explicitly: an empty/partial MD cohort (e.g. every
        # system failed equilibration) yields a labels frame without these
        # columns, and the bare label[[...]] selection would KeyError. A missing
        # MD result is a valid state -> just skip the metrics merge for it.
        merged = labels[["molecule_id", "md_candidate_id"]].merge(metrics, on="md_candidate_id", how="left")
        for row in merged.to_dict("records"):
            out.setdefault(row["molecule_id"], {}).update(row)
    if parameterization is not None:
        for row in parameterization.to_dict("records"):
            out.setdefault(row["molecule_id"], {})["parameterization_status"] = row.get("parameterization_status")
            out.setdefault(row["molecule_id"], {})["parameterization_backend"] = row.get("backend")
            out.setdefault(row["molecule_id"], {})["parameterization_warnings_json"] = row.get("warnings_json")
    return out


def _known_control_rows(inputs: dict[str, pd.DataFrame | None], detailed_count: int, summary_count: int) -> list[dict[str, Any]]:
    ranked = inputs["stage8_ranked"]
    aggregate = inputs["stage8_aggregate"]
    master = inputs["stage7_master"]
    md = _md_lookup(inputs)
    rows: list[dict[str, Any]] = []
    merged = ranked.merge(aggregate, on="molecule_id", how="left", suffixes=("_ranked", ""))
    if master is not None:
        keep = ["molecule_id", "inchi_key", "mw", "clogp", "tpsa", "hbd", "hba", "rotatable_bonds", "qed", "risk_flags_json"]
        merged = merged.merge(master[[c for c in keep if c in master.columns]], on="molecule_id", how="left", suffixes=("", "_master"))
    for idx, row in merged.sort_values("final_candidate_score", ascending=False).reset_index(drop=True).iterrows():
        md_row = md.get(row["molecule_id"], {})
        warnings = merge_warnings(_value(row, "warnings_json"), md_row.get("warnings_json"), md_row.get("parameterization_warnings_json"))
        label = "known_control_recovered" if row.get("screening_role") == "known_activity_reference" else "low_confidence_rejected"
        if _value(row, "best_pose_confidence", 0.0) < 0.25:
            label = "known_control_recovered" if row.get("screening_role") == "known_activity_reference" else "low_confidence_rejected"
        md_acceptance = md_row.get("md_acceptance_flag")
        md_label = str(md_row.get("md_stability_label", "")).lower()
        if md_acceptance is False or md_label == "md_unstable":
            label = "md_unstable_rejected"
            selection_reason = "Candidate completed Stage 11 MD but failed pose-stability thresholds; retained for audit display, not advancement."
        else:
            selection_reason = "Diagnostic known-control recovery row; retained for audit display, not a prospective prioritization claim."
        rows.append({
            "final_candidate_id": f"fcand_{idx + 1:03d}",
            "molecule_id": row["molecule_id"],
            "analog_id_if_available": None,
            "source": _value(row, "source", "unknown"),
            "subsource": _value(row, "subsource", _value(row, "subsource_ranked", "unknown")),
            "screening_role": _value(row, "screening_role", "unknown"),
            "standard_smiles": _value(row, "standard_smiles", _value(row, "standard_smiles_ranked")),
            "scaffold_id": _value(row, "scaffold_id", _value(row, "scaffold_id_ranked")),
            "novelty_bucket": _value(row, "novelty_bucket", _value(row, "novelty_bucket_ranked", "unknown")),
            "closest_known_molecule_id": _value(row, "closest_known_molecule_id", _value(row, "closest_known_molecule_id_ranked")),
            "tanimoto_to_closest_known": _value(row, "tanimoto_to_closest_known", _value(row, "tanimoto_to_closest_known_ranked")),
            "best_pose_id": _value(row, "best_screening_pose_id", _value(row, "best_pose_id")),
            "best_receptor_id": _value(row, "best_target_receptor_id"),
            "best_receptor_state": _value(row, "best_receptor_state"),
            "final_candidate_score": _value(row, "final_candidate_score", 0.0),
            "best_pose_confidence": _value(row, "best_pose_confidence", 0.0),
            "best_calibrated_confidence": _value(row, "best_calibrated_confidence", 0.0),
            "best_gnina_cnnscore": _value(row, "best_gnina_cnnscore"),
            "best_gnina_cnnaffinity": _value(row, "best_gnina_cnnaffinity"),
            "best_docking_score": _value(row, "best_docking_score"),
            "best_key_interaction_recall_consensus": _value(row, "best_key_interaction_recall_consensus"),
            "best_ifp_tanimoto_to_consensus": _value(row, "best_ifp_tanimoto_to_consensus"),
            "md_stability_label_if_available": md_row.get("md_stability_label", "not_available"),
            "md_key_interaction_persistence_if_available": md_row.get("interaction_persistence_component"),
            "medchem_risk_score": _value(row, "medchem_risk_score", 0.0),
            "risk_flags_json": _value(row, "risk_flags_json", _value(row, "medchem_flags_json", "[]")),
            "decision_label": label,
            "selected_for_detailed_dossier": idx < detailed_count,
            "selected_for_summary_table": idx < summary_count,
            "selection_reason": selection_reason,
            "nonclaim_statement": nonclaim_text(),
            "warnings_json": warnings,
        })
    return rows


def _rejected_analog_rows(inputs: dict[str, pd.DataFrame | None], start_index: int, summary_count: int) -> list[dict[str, Any]]:
    acceptance = inputs.get("stage9_analog_acceptance")
    analogs = inputs.get("stage9_analog_candidates")
    benchmark = inputs.get("stage10_analog_benchmark")
    if acceptance is None or acceptance.empty:
        return []
    rows: list[dict[str, Any]] = []
    df = acceptance[~acceptance["accepted_flag"].fillna(False)].copy()
    if analogs is not None:
        df = df.merge(analogs, on="analog_id", how="left", suffixes=("", "_analog"))
    if benchmark is not None:
        columns = ["analog_id", "scaffold_id", "novelty_bucket", "best_docking_score", "best_gnina_cnnscore", "best_gnina_cnnaffinity", "best_pose_confidence", "best_key_interaction_recall_consensus", "best_ifp_tanimoto_to_consensus", "medchem_risk_score", "warnings_json"]
        df = df.merge(benchmark[[c for c in columns if c in benchmark.columns]], on="analog_id", how="left", suffixes=("", "_bench"))
    for offset, row in df.head(max(summary_count - start_index, 0)).reset_index(drop=True).iterrows():
        idx = start_index + offset
        label = "score_hacking_rejected" if bool(_value(row, "score_hacking_flag", False)) else "low_confidence_rejected"
        rows.append({
            "final_candidate_id": f"fcand_{idx + 1:03d}",
            "molecule_id": _value(row, "parent_molecule_id", row["analog_id"]),
            "analog_id_if_available": row["analog_id"],
            "source": _value(row, "source", "stage9_generated_analog"),
            "subsource": _value(row, "strategy_name", "generated"),
            "screening_role": "generated_analog_negative_control",
            "standard_smiles": _value(row, "standard_smiles", _value(row, "analog_smiles")),
            "scaffold_id": _value(row, "scaffold_id"),
            "novelty_bucket": _value(row, "novelty_bucket", _value(row, "novelty_status", "analog")),
            "closest_known_molecule_id": _value(row, "closest_known_egfr_ligand", _value(row, "parent_molecule_id")),
            "tanimoto_to_closest_known": _value(row, "parent_tanimoto"),
            "best_pose_id": None,
            "best_receptor_id": None,
            "best_receptor_state": None,
            "final_candidate_score": _value(row, "analog_candidate_score", 0.0),
            "best_pose_confidence": _value(row, "analog_pose_confidence", _value(row, "best_pose_confidence", 0.0)),
            "best_calibrated_confidence": None,
            "best_gnina_cnnscore": _value(row, "analog_gnina_cnnscore", _value(row, "best_gnina_cnnscore")),
            "best_gnina_cnnaffinity": _value(row, "best_gnina_cnnaffinity"),
            "best_docking_score": _value(row, "best_docking_score"),
            "best_key_interaction_recall_consensus": _value(row, "analog_key_interaction_recall", _value(row, "best_key_interaction_recall_consensus")),
            "best_ifp_tanimoto_to_consensus": _value(row, "best_ifp_tanimoto_to_consensus"),
            "md_stability_label_if_available": "not_available",
            "md_key_interaction_persistence_if_available": None,
            "medchem_risk_score": _value(row, "medchem_risk_score", 0.0),
            "risk_flags_json": json.dumps([_value(row, "rejection_reason", "rejected_analog")]),
            "decision_label": label,
            "selected_for_detailed_dossier": False,
            "selected_for_summary_table": True,
            "selection_reason": f"Rejected generated analog retained as a negative-control example: {_value(row, 'rejection_reason', 'not accepted')}.",
            "nonclaim_statement": nonclaim_text(),
            "warnings_json": merge_warnings(_value(row, "warnings_json"), _value(row, "warnings_json_bench")),
        })
    return rows


def build_final_candidate_table(config_path: str | Path) -> dict[str, Any]:
    config = load_stage12_config(config_path)
    paths = stage12_paths(config)
    inputs = load_stage12_inputs(config)
    if inputs.get("stage8_ranked") is None or inputs.get("stage8_aggregate") is None:
        raise FileNotFoundError("Stage 12 requires Stage 8 ranked and aggregate candidate tables.")
    detailed_count = int(config["selection"]["detailed_dossier_count"])
    summary_count = int(config["selection"]["summary_table_count"])
    rows = _known_control_rows(inputs, detailed_count, summary_count)
    if config["selection"].get("include_rejected_analogs_for_negative_controls", True):
        rows.extend(_rejected_analog_rows(inputs, len(rows), summary_count))
    frame = pd.DataFrame(rows)
    for column in FINAL_CANDIDATE_COLUMNS:
        if column not in frame.columns:
            frame[column] = None
    frame = frame[FINAL_CANDIDATE_COLUMNS]
    validate_final_candidates(frame)
    ranked = frame.sort_values("final_candidate_score", ascending=False).reset_index(drop=True)
    ranked["final_rank_global"] = range(1, len(ranked) + 1)
    ranked["final_rank_within_source"] = ranked.groupby("source")["final_candidate_score"].rank(ascending=False, method="first").astype(int)
    ranked["final_rank_within_novelty_bucket"] = ranked.groupby("novelty_bucket")["final_candidate_score"].rank(ascending=False, method="first").astype(int)
    ranked["pose_confidence"] = ranked["best_pose_confidence"]
    ranked["calibrated_confidence"] = ranked["best_calibrated_confidence"]
    ranked["gnina_cnnscore"] = ranked["best_gnina_cnnscore"]
    ranked["cnnaffinity"] = ranked["best_gnina_cnnaffinity"]
    ranked["key_interaction_recall_consensus"] = ranked["best_key_interaction_recall_consensus"]
    ranked["ifp_tanimoto_to_consensus"] = ranked["best_ifp_tanimoto_to_consensus"]
    ranked["selected_for_dossier"] = ranked["selected_for_detailed_dossier"]
    ranked["recommended_next_action"] = ranked["decision_label"].map({
        "known_control_recovered": "manual_review_required",
        "score_hacking_rejected": "reject_score_hacking",
        "low_confidence_rejected": "manual_review_required",
        "md_unstable_rejected": "reject_md_unstable",
    }).fillna("experimental_testing_required")
    write_table(paths["processed"] / "final_candidate_selection.parquet", frame)
    write_table(paths["processed"] / "final_candidate_selection.csv", frame)
    ranked_out = ranked[FINAL_RANKED_COLUMNS]
    write_table(paths["processed"] / "final_ranked_candidates.parquet", ranked_out)
    write_table(paths["processed"] / "final_ranked_candidates.csv", ranked_out)
    table = paths["reports"] / "final_candidate_table"
    ranked_out.to_csv(table.with_suffix(".csv"), index=False)
    ranked_out.to_html(table.with_suffix(".html"), index=False, escape=True)
    try:
        ranked_out.to_excel(table.with_suffix(".xlsx"), index=False)
    except ImportError:
        pass
    return {"candidates": len(frame), "selected": int(frame["selected_for_detailed_dossier"].sum()), "table": str(paths["processed"] / "final_candidate_selection.parquet")}
