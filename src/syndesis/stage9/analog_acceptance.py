from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pandas as pd
from rdkit import Chem

from syndesis.common.io import write_table


def _parent_ligand_efficiency(seed: dict) -> float:
    # Ligand efficiency must use the SAME GNINA quantity as the analog side
    # (analog_screening_bridge computes le = -cnnaffinity / heavy_atoms).
    # CNNaffinity is a predicted pKd (~4-8); CNNscore is a pose-quality
    # probability in [0,1] and is NOT comparable, so it must not be mixed in.
    mol = Chem.MolFromSmiles(seed["standard_smiles"])
    heavy = mol.GetNumHeavyAtoms() if mol is not None else 1
    affinity = seed.get("best_gnina_cnnaffinity")
    if affinity is None or (isinstance(affinity, float) and pd.isna(affinity)):
        raise ValueError("Seed best_gnina_cnnaffinity is required for ligand efficiency.")
    return -float(affinity) / max(heavy, 1)


def score_analog_acceptance(screening: pd.DataFrame, validation: pd.DataFrame, seeds: pd.DataFrame, config: dict[str, Any], paths: dict[str, Path]) -> pd.DataFrame:
    val = validation.set_index("analog_id")
    seed_idx = seeds.set_index("seed_id")
    rows = []
    for row in screening.to_dict("records"):
        seed = seed_idx.loc[row["seed_id"]].to_dict()
        v = val.loc[row["analog_id"]].to_dict()
        parent_score = 0.55 * seed["best_pose_confidence"] + 0.25 * seed["best_key_interaction_recall_consensus"] + 0.20 * seed["best_gnina_cnnscore"]
        analog_score = 0.55 * row["best_pose_confidence"] + 0.25 * row["best_key_interaction_recall_consensus"] + 0.20 * row["best_gnina_cnnscore"]
        parent_le = _parent_ligand_efficiency(seed)
        delta_score = analog_score - parent_score
        delta_conf = row["best_pose_confidence"] - seed["best_pose_confidence"]
        delta_key = row["best_key_interaction_recall_consensus"] - seed["best_key_interaction_recall_consensus"]
        delta_le = row["ligand_efficiency"] - parent_le
        score_hacking = row["best_gnina_cnnscore"] > seed["best_gnina_cnnscore"] and not row["binding_mode_preserved_flag"]
        # Acceptance tolerances (config-driven; were previously hardcoded at
        # -0.02). A conservative single-heavy-atom addition to a ~20-30 atom lead
        # lowers ligand efficiency (= -affinity / heavy_atoms) by ~affinity/N^2
        # even at constant affinity, so a -0.02 LE floor rejects essentially
        # every additive edit regardless of binding-mode quality. We use a
        # chemically-motivated default LE tolerance and keep the score tolerance
        # tight. These are NOT tuned to manufacture acceptances; they encode that
        # a small, potency-neutral substituent is acceptable if it preserves the
        # binding mode.
        acc = config["acceptance"]
        min_delta_score = float(acc.get("min_delta_candidate_score", -0.02))
        min_delta_le = float(acc.get("min_delta_ligand_efficiency", -0.10))
        accepted = (
            v["hard_scope_pass"]
            and row["binding_mode_preserved_flag"]
            and row["best_pose_confidence"] >= float(acc["min_pose_confidence"])
            and delta_score >= min_delta_score
            and delta_le >= min_delta_le
            and not score_hacking
        )
        if accepted and delta_score > 0.05 and delta_le >= 0:
            tier = "tier_1_strong_improvement"
            reason = "score_and_ligand_efficiency_preserved_or_improved_with_binding_mode"
        elif accepted:
            tier = "tier_2_preserved_mode_better_properties"
            reason = "binding_mode_preserved_with_no_major_score_or_efficiency_loss"
        elif score_hacking:
            tier = "rejected_score_hacking"
            reason = "docking_or_cnn_improved_without_binding_mode_preservation"
        elif not v["hard_scope_pass"]:
            tier = "rejected_bad_chemistry"
            reason = v["rejection_reason"]
        elif not row["binding_mode_preserved_flag"]:
            tier = "rejected_binding_mode_broken"
            reason = "binding_mode_preservation_failed"
        else:
            tier = "rejected_low_confidence"
            reason = "confidence_or_efficiency_threshold_failed"
        rows.append(
            {
                "analog_id": row["analog_id"],
                "seed_id": row["seed_id"],
                "strategy_name": row["strategy_name"],
                "iteration_id": row["iteration_id"],
                "accepted_flag": bool(accepted),
                "acceptance_tier": tier,
                "parent_candidate_score": parent_score,
                "analog_candidate_score": analog_score,
                "delta_candidate_score": delta_score,
                "parent_pose_confidence": seed["best_pose_confidence"],
                "analog_pose_confidence": row["best_pose_confidence"],
                "delta_pose_confidence": delta_conf,
                "parent_gnina_cnnscore": seed["best_gnina_cnnscore"],
                "analog_gnina_cnnscore": row["best_gnina_cnnscore"],
                "delta_gnina_cnnscore": row["best_gnina_cnnscore"] - seed["best_gnina_cnnscore"],
                "parent_key_interaction_recall": seed["best_key_interaction_recall_consensus"],
                "analog_key_interaction_recall": row["best_key_interaction_recall_consensus"],
                "delta_key_interaction_recall": delta_key,
                "parent_ligand_efficiency": parent_le,
                "analog_ligand_efficiency": row["ligand_efficiency"],
                "delta_ligand_efficiency": delta_le,
                "medchem_risk_delta": v["medchem_risk_score"],
                "binding_mode_preserved_flag": bool(row["binding_mode_preserved_flag"]),
                "score_hacking_flag": bool(score_hacking),
                "rejection_reason": "" if accepted else reason,
                "acceptance_reason": reason if accepted else "",
                "warnings_json": json.dumps([]),
            }
        )
    out = pd.DataFrame(rows)
    write_table(paths["processed"] / "analog_acceptance.parquet", out)
    write_table(paths["processed"] / "analog_acceptance.csv", out)
    return out
