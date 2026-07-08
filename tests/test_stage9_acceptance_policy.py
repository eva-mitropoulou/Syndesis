import pandas as pd

from egfr_dockingforge.stage9.analog_acceptance import score_analog_acceptance


def test_docking_score_improvement_alone_is_insufficient(tmp_path):
    seeds = pd.DataFrame(
        [
            {
                "seed_id": "seed_001",
                "standard_smiles": "Cc1ccccc1",
                "best_pose_confidence": 0.6,
                "best_key_interaction_recall_consensus": 0.5,
                "best_gnina_cnnscore": 0.4,
            }
        ]
    )
    validation = pd.DataFrame([{"analog_id": "a1", "hard_scope_pass": True, "rejection_reason": "", "medchem_risk_score": 0.0}])
    screening = pd.DataFrame(
        [
            {
                "analog_id": "a1",
                "seed_id": "seed_001",
                "strategy_name": "rdkit_rule_based",
                "iteration_id": "iter_001",
                "best_pose_confidence": 0.7,
                "best_key_interaction_recall_consensus": 0.0,
                "best_gnina_cnnscore": 0.9,
                "ligand_efficiency": 1.0,
                "binding_mode_preserved_flag": False,
            }
        ]
    )
    config = {"acceptance": {"min_pose_confidence": 0.3}}
    out = score_analog_acceptance(screening, validation, seeds, config, {"processed": tmp_path})
    assert not bool(out["accepted_flag"].item())
    assert out["acceptance_tier"].item() == "rejected_score_hacking"


def test_identical_to_parent_analog_is_accepted(tmp_path):
    """Guards against the ligand-efficiency unit-mismatch bug: an analog that is
    identical to its parent (same scores, binding mode preserved) must pass. The
    parent LE (from cnnaffinity) and analog LE (from cnnaffinity) must be on the
    same scale, so delta_le ~= 0 rather than a large spurious negative."""
    seeds = pd.DataFrame(
        [
            {
                "seed_id": "seed_001",
                "standard_smiles": "Cc1ccccc1",
                "best_pose_confidence": 0.7,
                "best_key_interaction_recall_consensus": 0.6,
                "best_gnina_cnnscore": 0.8,
                "best_gnina_cnnaffinity": 6.0,
            }
        ]
    )
    validation = pd.DataFrame([{"analog_id": "a1", "hard_scope_pass": True, "rejection_reason": "", "medchem_risk_score": 0.0}])
    # heavy atom count of toluene (Cc1ccccc1) is 7; parent le = -6.0/7 = -0.857.
    parent_le = -6.0 / 7
    screening = pd.DataFrame(
        [
            {
                "analog_id": "a1",
                "seed_id": "seed_001",
                "strategy_name": "rdkit_rule_based",
                "iteration_id": "iter_001",
                "best_pose_confidence": 0.7,
                "best_key_interaction_recall_consensus": 0.6,
                "best_gnina_cnnscore": 0.8,
                "ligand_efficiency": parent_le,
                "binding_mode_preserved_flag": True,
            }
        ]
    )
    config = {"acceptance": {"min_pose_confidence": 0.3}}
    out = score_analog_acceptance(screening, validation, seeds, config, {"processed": tmp_path})
    assert abs(float(out["delta_ligand_efficiency"].item())) < 1e-6
    assert bool(out["accepted_flag"].item())


def test_parent_ligand_efficiency_uses_cnnaffinity_not_cnnscore():
    """The parent LE numerator must be cnnaffinity (pKd scale), not cnnscore."""
    from egfr_dockingforge.stage9.analog_acceptance import _parent_ligand_efficiency

    seed = {"standard_smiles": "Cc1ccccc1", "best_gnina_cnnscore": 0.8, "best_gnina_cnnaffinity": 6.0}
    le = _parent_ligand_efficiency(seed)
    assert abs(le - (-6.0 / 7)) < 1e-6
