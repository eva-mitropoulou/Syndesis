from __future__ import annotations

STAGE9_TO_STAGE10 = {
    "rdkit_rule_based": "rdkit_rule_based",
    "reinvent4_baseline": "reinvent4_baseline_if_available",
    "single_agent": "single_agent",
    "council_no_feedback": "council_loop",
    "council_with_tool_feedback": "council_plus_prolif",
    "council_interaction_constrained": "council_plus_prolif_pose_confidence",
}

STRATEGIES = [
    "random_analog_enumeration",
    "docking_score_only_optimization",
    "gnina_only_optimization",
    "rdkit_rule_based",
    "reinvent4_baseline_if_available",
    "single_agent",
    "council_loop",
    "council_plus_prolif",
    "council_plus_prolif_pose_confidence",
    "council_plus_prolif_pose_confidence_md",
]


def normalize_strategy_name(name: str) -> str:
    return STAGE9_TO_STAGE10.get(name, name)
