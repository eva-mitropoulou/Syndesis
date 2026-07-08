import pandas as pd

from egfr_dockingforge.stage10.benchmark_manifest import build_budget_audit


def test_budget_violation_recorded_for_under_generation(tmp_path):
    master = pd.DataFrame(
        [
            {
                "strategy_id": "strat_01",
                "seed_id": "seed_001",
                "valid_molecule_flag": True,
                "unique_flag": True,
                "best_pose_confidence": 0.4,
            }
        ]
    )
    manifest = pd.DataFrame(
        [
            {
                "strategy_id": "strat_01",
                "generation_budget": 6,
                "screening_budget": 3,
                "enabled_flag": True,
            }
        ]
    )
    inputs = {"stage9_seeds": pd.DataFrame([{"seed_id": "seed_001"}])}
    out = build_budget_audit(master, manifest, inputs, {}, {"processed": tmp_path})
    assert out["budget_violation_reason"].item() == "generated_fewer_than_budget"
