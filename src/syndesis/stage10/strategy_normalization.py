from __future__ import annotations

STAGE9_TO_STAGE10 = {
    "rdkit_rule_based": "rdkit_rule_based",
}

STRATEGIES = [
    "random_analog_enumeration",
    "docking_score_only_optimization",
    "gnina_only_optimization",
    "rdkit_rule_based",
    "rdkit_rule_based_md",
]


def normalize_strategy_name(name: str) -> str:
    return STAGE9_TO_STAGE10.get(name, name)
