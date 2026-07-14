from syndesis.stage10.schemas import TABLE_SCHEMAS


def test_stage10_schemas_have_required_tables():
    required = {
        "ablation_strategy_manifest",
        "strategy_budget_audit",
        "analog_benchmark_master",
        "score_hacking_cases",
        "strategy_metrics",
        "seed_strategy_metrics",
        "statistical_comparisons",
        "ablation_summary",
        "diversity_novelty_metrics",
        "compute_cost_metrics",
    }
    assert required.issubset(TABLE_SCHEMAS)
    assert all(len(cols) == len(set(cols)) for cols in TABLE_SCHEMAS.values())
