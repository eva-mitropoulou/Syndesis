from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pandas as pd

from syndesis.common.io import write_table


def _row(base: dict[str, Any], receptor: str, pose: str, ligand_file: str) -> dict[str, Any]:
    """Assemble a full MD-candidate-manifest row from a partial candidate dict."""
    return {
        "molecule_id": base["molecule_id"],
        "source": base.get("source", ""),
        "parent_molecule_id": base.get("parent_molecule_id", base["molecule_id"]),
        "candidate_origin_stage": base["candidate_origin_stage"],
        "standard_smiles": base.get("standard_smiles", ""),
        "prepared_smiles": base.get("standard_smiles", ""),
        "best_pose_id": pose,
        "receptor_id": receptor,
        "receptor_state": base.get("receptor_state", ""),
        "protein_ligand_complex_file": "",
        "ligand_file": ligand_file,
        "pose_file": base.get("pose_file", ""),
        "receptor_file": f"data/processed/stage2/ensemble_receptors/{receptor}.pdb",
        "stage8_candidate_score": base.get("stage8_candidate_score", 0.0),
        "stage9_acceptance_tier": base.get("stage9_acceptance_tier", ""),
        "stage10_strategy_name": base.get("stage10_strategy_name", ""),
        "pose_confidence": base.get("pose_confidence", 0.0),
        "gnina_cnnscore": base.get("gnina_cnnscore", 0.0),
        "key_interaction_recall_consensus": base.get("key_interaction_recall_consensus", 0.0),
        "novelty_bucket": base.get("novelty_bucket", ""),
        "selection_reason": base.get("selection_reason", ""),
    }


def _accepted_analog_candidates(inputs: dict[str, pd.DataFrame | None]) -> list[dict[str, Any]]:
    """Accepted Stage-9 analogs eligible for MD, joined with their screening pose."""
    acceptance = inputs.get("stage9_acceptance")
    screening = inputs.get("stage9_screening")
    if acceptance is None or acceptance.empty or "accepted_flag" not in acceptance.columns:
        return []
    accepted = acceptance[acceptance["accepted_flag"].fillna(False).astype(bool)]
    if accepted.empty:
        return []
    screen_idx = screening.set_index("analog_id") if screening is not None and not screening.empty else None
    out = []
    for rec in accepted.to_dict("records"):
        analog_id = rec.get("analog_id")
        scr: dict[str, Any] = {}
        if screen_idx is not None and analog_id in screen_idx.index:
            s = screen_idx.loc[analog_id]
            scr = (s.iloc[0] if isinstance(s, pd.DataFrame) else s).to_dict()
        pose = scr.get("best_pose_id", "")
        out.append(
            {
                "molecule_id": analog_id,
                "source": f"stage9_analog_{rec.get('strategy_name', '')}",
                "parent_molecule_id": rec.get("seed_id", ""),
                "candidate_origin_stage": "stage9_accepted_analog",
                "standard_smiles": rec.get("analog_smiles") or scr.get("standard_smiles", ""),
                "receptor_id": scr.get("best_receptor_id", ""),
                "receptor_state": scr.get("best_receptor_state", ""),
                "ligand_file": f"data/processed/stage9/prepared_analogs/{analog_id}.sdf",
                "pose_file": f"data/processed/stage9/stage8_mini_screen/gnina_ligands/{pose}.pdb",
                "best_pose_id": pose,
                "stage8_candidate_score": scr.get("best_pose_confidence", 0.0),
                "stage9_acceptance_tier": rec.get("acceptance_tier", ""),
                "stage10_strategy_name": rec.get("strategy_name", ""),
                "pose_confidence": scr.get("best_pose_confidence", 0.0),
                "gnina_cnnscore": scr.get("best_gnina_cnnscore", 0.0),
                "key_interaction_recall_consensus": scr.get("best_key_interaction_recall_consensus", 0.0),
                "novelty_bucket": "stage9_accepted_analog",
                "selection_reason": f"stage9_accepted_analog:{rec.get('acceptance_tier','')}",
            }
        )
    return out


def select_md_candidates(inputs: dict[str, pd.DataFrame | None], config: dict[str, Any], paths: dict[str, Path]) -> pd.DataFrame:
    ranked = inputs["stage8_ranked_candidates"].copy()
    max_n = int(config["selection"]["max_quick_md_candidates"])
    include_controls = bool(config["selection"].get("include_known_controls", True))
    # Finalist policy: known controls are always finalists (diagnostic anchors);
    # the top-K accepted analogs (by pose confidence) are promoted to finalist and
    # receive the full replicated, hardened MD. The remaining accepted analogs are
    # retained in the manifest as "selected_md_pending" (reported as future work)
    # so the scope is explicit and honest.
    max_finalist_analogs = int(config["selection"].get("max_finalist_analog_candidates", 3))
    replicates = int(config.get("md", {}).get("finalist_replicates", 3))

    controls: list[dict[str, Any]] = []
    if include_controls:
        for row in ranked.head(max_n).to_dict("records"):
            controls.append(
                {
                    "molecule_id": row["molecule_id"],
                    "source": row.get("source", ""),
                    "parent_molecule_id": row.get("closest_known_molecule_id", row["molecule_id"]),
                    "candidate_origin_stage": "stage8_known_control_or_candidate",
                    "standard_smiles": row["standard_smiles"],
                    "receptor_id": row.get("best_target_receptor_id", ""),
                    "receptor_state": row.get("best_receptor_state", ""),
                    "ligand_file": f"data/processed/stage7/prepared_ligands/prep_{row['molecule_id']}.sdf",
                    "pose_file": f"data/processed/stage8/gnina_ligands/{row.get('best_screening_pose_id', '')}.pdb",
                    "best_pose_id": row.get("best_screening_pose_id", ""),
                    "stage8_candidate_score": row.get("final_candidate_score", 0.0),
                    "pose_confidence": row.get("best_pose_confidence", 0.0),
                    "gnina_cnnscore": row.get("best_gnina_cnnscore", 0.0),
                    "key_interaction_recall_consensus": row.get("best_key_interaction_recall_consensus", 0.0),
                    "novelty_bucket": row.get("novelty_bucket", ""),
                    "selection_reason": "known_control_or_top_stage8_candidate",
                    "_is_finalist": True,
                }
            )

    analogs: list[dict[str, Any]] = []
    if config["selection"].get("require_stage9_accepted_for_analogs", True):
        analogs = _accepted_analog_candidates(inputs)
        # Rank accepted analogs by pose confidence (then GNINA CNNscore) and
        # promote the top-K to finalist.
        analogs.sort(key=lambda c: (c.get("pose_confidence", 0.0), c.get("gnina_cnnscore", 0.0)), reverse=True)
        for j, cand in enumerate(analogs):
            finalist = j < max_finalist_analogs
            cand["_is_finalist"] = finalist
            cand["selection_reason"] = (
                f"finalist_top_accepted_analog:{cand.get('stage9_acceptance_tier','')}" if finalist
                else f"selected_md_pending_future_work:{cand.get('stage9_acceptance_tier','')}"
            )

    rows = []
    for i, cand in enumerate([*controls, *analogs], start=1):
        base = _row(cand, cand.get("receptor_id", ""), cand.get("best_pose_id", ""), cand.get("ligand_file", ""))
        base["md_candidate_id"] = f"mdcand_{i:03d}"
        finalist = bool(cand.get("_is_finalist", False))
        base["is_md_finalist"] = finalist
        base["selected_for_quick_md"] = finalist
        base["selected_for_replicate_md"] = finalist
        base["planned_replicates"] = replicates if finalist else 0
        base["md_status"] = "finalist_scheduled" if finalist else "selected_md_pending"
        base["warnings_json"] = json.dumps([] if finalist else ["selected_for_md_pending_future_work"])
        rows.append(base)
    out = pd.DataFrame(rows)
    write_table(paths["processed"] / "md_candidate_manifest.parquet", out)
    write_table(paths["processed"] / "md_candidate_manifest.csv", out)
    return out
