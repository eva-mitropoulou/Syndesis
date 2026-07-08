import json

from egfr_dockingforge.stage9.agent_schemas import REQUIRED_PROPOSAL_FIELDS, parse_agent_json


def test_valid_agent_output_parses_without_chain_of_thought():
    payload = {field: "" for field in REQUIRED_PROPOSAL_FIELDS}
    payload.update({"requires_tool": True, "confidence": 0.7, "rationale_summary": "Concise tool-verifiable rationale."})
    result = parse_agent_json(json.dumps(payload))
    assert result.schema_valid
    assert result.proposal_status == "schema_valid"


def test_invalid_agent_output_fails_clearly():
    result = parse_agent_json("not json")
    assert not result.schema_valid
    assert result.proposal_status == "failed_agent_output"


def test_chain_of_thought_field_rejected():
    payload = {field: "" for field in REQUIRED_PROPOSAL_FIELDS}
    payload["chain_of_thought"] = "hidden reasoning"
    result = parse_agent_json(json.dumps(payload))
    assert not result.schema_valid
