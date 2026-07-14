from __future__ import annotations


def accepted_pre_md(
    row: dict,
    min_pose_confidence: float = 0.30,
    max_medchem_risk: float = 0.40,
    min_delta_candidate_score: float = -0.02,
    min_delta_ligand_efficiency: float = -0.10,
) -> bool:
    # Tolerances must match Stage 9 analog_acceptance (config-driven) so the
    # Stage 10 ablation reflects the same acceptances rather than re-deriving
    # them with a stale hardcoded gate.
    return bool(
        row.get("valid_molecule_flag")
        and row.get("unique_flag")
        and row.get("hard_scope_pass")
        and not row.get("covalent_warhead_flag")
        and not row.get("reactive_flag")
        and row.get("medchem_risk_score", 1.0) <= max_medchem_risk
        and row.get("binding_mode_preserved_flag")
        and row.get("best_pose_confidence", 0.0) >= min_pose_confidence
        and row.get("delta_candidate_score", -999.0) >= min_delta_candidate_score
        and row.get("delta_ligand_efficiency", -999.0) >= min_delta_ligand_efficiency
    )


def score_hacking(row: dict) -> bool:
    improved_cnn = row.get("delta_gnina_cnnscore", 0.0) > 0
    improved_score = row.get("delta_candidate_score", 0.0) > 0
    worsened_pose = row.get("delta_pose_confidence", 0.0) < 0
    worsened_key = row.get("delta_key_interaction_recall", 0.0) < 0
    broken_mode = not bool(row.get("binding_mode_preserved_flag"))
    bad_risk = row.get("medchem_risk_score", 0.0) > 0.40
    return bool((improved_cnn and (worsened_pose or worsened_key or broken_mode)) or (improved_score and bad_risk))
