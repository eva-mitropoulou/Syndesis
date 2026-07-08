from __future__ import annotations

import json

AGENT_ROLES = {
    "medicinal_chemist": "Propose conservative, synthetically plausible transformations while preserving the Stage 8 scaffold.",
    "docking_scientist": "Use receptor-state and pose evidence to avoid transformations that are likely to break the docked binding mode.",
    "interaction_analyst": "Protect key ProLIF interactions and propose recoverable interaction improvements only.",
    "admet_safety_critic": "Flag property, PAINS, reactivity, and lipophilicity risks without making clinical claims.",
    "skeptical_reviewer": "Reject score-hacking and proposals unsupported by tool-verifiable evidence.",
    "orchestrator": "Consolidate transformation requests and never invent final molecules directly.",
}


def compact_prompt(role: str, seed_summary: dict, edit_sites: list[dict]) -> str:
    seed_id = str(seed_summary.get("seed_id", "seed_unknown"))
    edit_site_id = str(edit_sites[0].get("edit_site_id", f"{seed_id}_site_01")) if edit_sites else f"{seed_id}_site_01"
    allowed_classes = [
        "small_substituent_scan",
        "halogen_scan",
        "heteroatom_swap",
        "solubilizing_tail_tuning",
        "conservative_bioisostere",
        "brics_recombination",
    ]
    skeleton = {
        "proposal_id": f"proposal_{seed_id}_{role}",
        "seed_id": seed_id,
        "agent_role": role,
        "proposed_transformation_class": "small_substituent_scan",
        "edit_site_id": edit_site_id,
        "transformation_description": "one tool-executable analog transformation",
        "expected_effect": "short expected binding or developability effect",
        "expected_preserved_interactions": [],
        "expected_new_interactions": [],
        "medchem_risk_prediction": "low",
        "priority": 1,
        "requires_tool": True,
        "tool_name": "rdkit_transform",
        "tool_arguments_json": {},
        "rationale_summary": "concise tool-verifiable rationale",
        "confidence": 0.5,
        "reject_if_conditions_json": ["invalid_valence", "binding_mode_broken"],
    }
    return (
        "Return exactly one compact JSON object. Do not output chain-of-thought, markdown, wrapper keys, "
        "or nested top-level objects. Use every top-level key from the JSON skeleton exactly once. "
        f"Role: {role}. Task: propose one conservative, tool-executable analog transformation. "
        f"Allowed transformation classes: {allowed_classes}. "
        f"Seed: {seed_summary}. Editable sites: {edit_sites}. "
        f"Fill this exact flat JSON shape with project-specific values: {json.dumps(skeleton, sort_keys=True)}"
    )
