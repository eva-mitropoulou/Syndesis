import pandas as pd

from egfr_dockingforge.stage9.baseline_benchmark import benchmark_strategies


def test_benchmark_computes_accepted_rate(tmp_path):
    candidates = pd.DataFrame([{"analog_id": "a1", "strategy_name": "rdkit_rule_based", "standard_smiles": "C", "seed_id": "s"}])
    validation = pd.DataFrame([{"analog_id": "a1", "valid_molecule_flag": True, "hard_scope_pass": True}])
    screening = pd.DataFrame([{"analog_id": "a1", "strategy_name": "rdkit_rule_based"}])
    acceptance = pd.DataFrame(
        [
            {
                "analog_id": "a1",
                "strategy_name": "rdkit_rule_based",
                "accepted_flag": True,
                "delta_candidate_score": 0.1,
                "delta_pose_confidence": 0.1,
                "delta_key_interaction_recall": 0.1,
                "score_hacking_flag": False,
                "binding_mode_preserved_flag": True,
            }
        ]
    )
    config = {"loop": {"strategies": ["rdkit_rule_based"], "max_iterations": 1}}
    out = benchmark_strategies(candidates, validation, screening, acceptance, config, {"processed": tmp_path})
    assert out["accepted_analog_rate"].item() == 1.0
