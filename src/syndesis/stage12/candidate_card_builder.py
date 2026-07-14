from __future__ import annotations

import hashlib
import json
import subprocess
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import pandas as pd

from syndesis.common.io import write_json
from syndesis.stage12.candidate_card_schema import validate_candidate_card
from syndesis.stage12.evidence_aggregation import best_pose_file, load_selection_and_inputs
from syndesis.stage12.nonclaim_generator import nonclaim_dict
from syndesis.stage12.risk_summary import main_risks, parse_json_list
from syndesis.stage12.structure_figures import write_candidate_structure_svg


def _git_commit() -> str | None:
    try:
        return subprocess.check_output(["git", "rev-parse", "HEAD"], text=True).strip()
    except Exception:
        return None


def _hash_file(path: str | Path) -> str | None:
    p = Path(path)
    if not p.exists():
        return None
    return hashlib.sha256(p.read_bytes()).hexdigest()


def _master_lookup(inputs: dict[str, pd.DataFrame | None]) -> dict[str, dict[str, Any]]:
    master = inputs.get("stage7_master")
    if master is None:
        return {}
    return {row["molecule_id"]: row for row in master.to_dict("records")}


def _analog_lookup(inputs: dict[str, pd.DataFrame | None]) -> dict[str, dict[str, Any]]:
    analogs = inputs.get("stage9_analog_candidates")
    if analogs is None:
        return {}
    return {row["analog_id"]: row for row in analogs.to_dict("records")}


def card_from_row(row: dict[str, Any], inputs: dict[str, pd.DataFrame | None], config_path: str | Path) -> dict[str, Any]:
    master = _master_lookup(inputs).get(row["molecule_id"], {})
    analog = _analog_lookup(inputs).get(row.get("analog_id_if_available"), {})
    risks = main_risks(row)
    pose_file = best_pose_file(row.get("best_pose_id"))
    card = {
        "card_version": "stage12.v1",
        "final_candidate_id": row["final_candidate_id"],
        "molecule_id": row["molecule_id"],
        "source": row["source"],
        "subsource": row["subsource"],
        "screening_role": row["screening_role"],
        "standard_smiles": row["standard_smiles"],
        "inchi_key": master.get("inchi_key") or analog.get("inchi_key"),
        "scaffold_id": row["scaffold_id"],
        "novelty_bucket": row["novelty_bucket"],
        "closest_known_egfr_ligand": {
            "molecule_id": row.get("closest_known_molecule_id"),
            "name_if_available": None,
            "similarity": row.get("tanimoto_to_closest_known"),
            "known_activity_summary_if_available": master.get("known_activity_status"),
            "source": "stage7_candidate_library",
        },
        "parent_analog_lineage": {
            "parent_molecule_id": analog.get("parent_molecule_id"),
            "parent_smiles": analog.get("parent_smiles"),
            "transformation_class": analog.get("transformation_class"),
            "seed_id": analog.get("seed_id"),
            "strategy_name": analog.get("strategy_name"),
        },
        "best_pose": {
            "pose_id": row.get("best_pose_id"),
            "receptor_id": row.get("best_receptor_id"),
            "receptor_state": row.get("best_receptor_state"),
            "pose_file": pose_file,
            "receptor_file": f"data/processed/stage8/prolif_receptors/ensemble_receptors__{row.get('best_receptor_id')}.h.pdb" if row.get("best_receptor_id") else None,
            "ligand_file": pose_file,
            "missing_pose_reason": None if pose_file else "pose file unavailable or candidate is a rejected analog example",
        },
        "scores": {
            "docking_score": row.get("best_docking_score"),
            "gnina_cnnscore": row.get("best_gnina_cnnscore"),
            "gnina_cnnaffinity": row.get("best_gnina_cnnaffinity"),
            "pose_confidence": row.get("best_pose_confidence"),
            "calibrated_confidence": row.get("best_calibrated_confidence"),
            "final_candidate_score": row.get("final_candidate_score"),
        },
        "interactions": {
            "key_interactions_recovered": [],
            "key_interactions_missing": [],
            "ifp_tanimoto_to_consensus": row.get("best_ifp_tanimoto_to_consensus"),
            "key_interaction_recall_consensus": row.get("best_key_interaction_recall_consensus"),
            "prolif_summary": "Interaction evidence is summarized from Stage 8 ProLIF-derived tables.",
        },
        "md": {
            "status": "not_available" if row.get("md_stability_label_if_available") in (None, "not_available", "md_failed_setup") else "available",
            "ligand_rmsd_summary": None,
            "key_interaction_persistence": row.get("md_key_interaction_persistence_if_available"),
            "pose_stability_label": row.get("md_stability_label_if_available"),
            "replicate_count": 0 if row.get("md_stability_label_if_available") in (None, "not_available", "md_failed_setup") else 1,
        },
        "medchem": {
            "mw": master.get("mw"),
            "clogp": master.get("clogp"),
            "tpsa": master.get("tpsa"),
            "hbd": master.get("hbd"),
            "hba": master.get("hba"),
            "rotatable_bonds": master.get("rotatable_bonds"),
            "qed": master.get("qed"),
            "risk_flags": parse_json_list(row.get("risk_flags_json")),
        },
        "evidence_summary": {
            "selected_because": row.get("selection_reason"),
            "main_risks": risks,
            "rejection_or_caution_reason": "; ".join(risks) if risks else "no additional risk flags recorded",
            "what_would_need_experimental_testing": "Biochemical potency, cellular response, selectivity, PK, toxicity, and structural confirmation.",
        },
        "non_claims": nonclaim_dict(),
        "provenance": {
            "input_tables": list(inputs.keys()),
            "model_versions": {"pose_confidence": "stage6_selected_model_if_available"},
            "tool_versions": {"rdkit": "used_if_installed_for_2d_depictions"},
            "config_hashes": {str(config_path): _hash_file(config_path)},
            "generation_timestamp": datetime.now(timezone.utc).isoformat(),
            "git_commit_if_available": _git_commit(),
            "environment_file": "environment.yml",
        },
    }
    validate_candidate_card(card)
    return card


def build_candidate_cards(config_path: str | Path) -> dict[str, Any]:
    _, paths, inputs, selection = load_selection_and_inputs(config_path)
    count = 0
    for row in selection.to_dict("records"):
        card = card_from_row(row, inputs, config_path)
        write_json(paths["cards"] / f"{row['final_candidate_id']}.json", card)
        write_candidate_structure_svg(row, paths["figures"])
        count += 1
    return {"cards": count, "cards_dir": str(paths["cards"])}
