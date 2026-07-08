from __future__ import annotations

import json
from dataclasses import dataclass
from typing import Any

REQUIRED_PROPOSAL_FIELDS = [
    "proposal_id",
    "seed_id",
    "agent_role",
    "proposed_transformation_class",
    "edit_site_id",
    "transformation_description",
    "expected_effect",
    "expected_preserved_interactions",
    "expected_new_interactions",
    "medchem_risk_prediction",
    "priority",
    "requires_tool",
    "tool_name",
    "tool_arguments_json",
    "rationale_summary",
    "confidence",
    "reject_if_conditions_json",
]


@dataclass(frozen=True)
class AgentParseResult:
    parsed_json: dict[str, Any] | None
    schema_valid: bool
    repair_attempts: int
    proposal_status: str
    warning: str


def _strip_json_fences(text: str) -> str:
    """Remove markdown code fences and locate the outermost JSON object.

    vLLM guided-JSON (response_format=json_object) returns bare JSON, but some
    providers/templates wrap it in ```json ... ``` fences or add leading prose.
    Stripping these prevents a valid proposal from being discarded on a parse
    error.
    """
    stripped = text.strip()
    if stripped.startswith("```"):
        # Drop the opening fence line (``` or ```json) and any trailing fence.
        lines = stripped.splitlines()
        if lines and lines[0].startswith("```"):
            lines = lines[1:]
        if lines and lines[-1].strip().startswith("```"):
            lines = lines[:-1]
        stripped = "\n".join(lines).strip()
    # Fall back to the substring between the first '{' and last '}'.
    if not stripped.startswith("{"):
        start = stripped.find("{")
        end = stripped.rfind("}")
        if start != -1 and end != -1 and end > start:
            stripped = stripped[start : end + 1]
    return stripped


def parse_agent_json(raw_text: str, repair_attempts: int = 0) -> AgentParseResult:
    try:
        payload = json.loads(raw_text)
    except json.JSONDecodeError:
        try:
            payload = json.loads(_strip_json_fences(raw_text))
        except json.JSONDecodeError as exc:
            return AgentParseResult(None, False, repair_attempts, "failed_agent_output", str(exc))
    missing = [field for field in REQUIRED_PROPOSAL_FIELDS if field not in payload]
    forbidden = [field for field in payload if "chain_of_thought" in field.lower()]
    if missing or forbidden:
        return AgentParseResult(
            payload,
            False,
            repair_attempts,
            "failed_agent_output",
            f"missing={missing}; forbidden={forbidden}",
        )
    return AgentParseResult(payload, True, repair_attempts, "schema_valid", "")


def proposal_json_schema() -> dict[str, Any]:
    return {
        "type": "object",
        "required": REQUIRED_PROPOSAL_FIELDS,
        "additionalProperties": True,
        "properties": {field: {"type": ["string", "number", "boolean", "array", "object"]} for field in REQUIRED_PROPOSAL_FIELDS},
    }
