from __future__ import annotations

import json
import shutil
from dataclasses import dataclass
from typing import Any

import requests


@dataclass(frozen=True)
class LLMResult:
    ok: bool
    raw_text: str
    status: str
    warnings: list[str]
    token_count: int = 0


def provider_available(config: dict[str, Any]) -> tuple[bool, str]:
    provider = config["llm"]["provider"]
    if provider == "disabled_mock_for_tests":
        return True, "disabled_mock_for_tests"
    if provider == "ollama":
        if not shutil.which("ollama"):
            return False, "ollama_cli_not_on_path"
        try:
            response = requests.get(f"{config['llm']['base_url'].rstrip('/')}/api/tags", timeout=5)
        except requests.RequestException as exc:
            return False, f"ollama_api_unreachable:{exc.__class__.__name__}"
        return response.ok, f"ollama_api_status_{response.status_code}"
    if provider in {"llama_cpp_server", "vllm", "openai_compatible_local"}:
        try:
            response = requests.get(
                config["llm"]["base_url"].rstrip("/") + "/v1/models",
                headers=_headers(config),
                timeout=5,
            )
        except requests.RequestException as exc:
            return False, f"{provider}_api_unreachable:{exc.__class__.__name__}"
        return response.ok, f"{provider}_api_status_{response.status_code}"
    return False, f"unsupported_provider:{provider}"


def _headers(config: dict[str, Any]) -> dict[str, str]:
    api_key = config.get("llm", {}).get("api_key")
    if not api_key:
        return {}
    return {"Authorization": f"Bearer {api_key}"}


def _json_payload(prompt: str, config: dict[str, Any]) -> dict[str, Any]:
    llm = config["llm"]
    payload = {
        "model": llm["model_name"],
        "messages": [
            {
                "role": "system",
                "content": (
                    "Return exactly one compact JSON object matching the requested flat schema. "
                    "Do not include markdown, wrapper keys, nested top-level objects, or chain-of-thought."
                ),
            },
            {"role": "user", "content": prompt},
        ],
        "temperature": float(llm.get("temperature", 0.2)),
        "top_p": float(llm.get("top_p", 0.9)),
        "max_tokens": int(llm.get("max_tokens", 2048)),
        "response_format": {"type": "json_object"},
    }
    chat_template_kwargs = llm.get("chat_template_kwargs")
    if chat_template_kwargs:
        payload["chat_template_kwargs"] = chat_template_kwargs
    return payload


def _extract_content(payload: dict[str, Any]) -> tuple[str, int]:
    choices = payload.get("choices") or []
    if not choices:
        return "", int((payload.get("usage") or {}).get("total_tokens") or 0)
    message = choices[0].get("message") or {}
    content = message.get("content") or ""
    return str(content).strip(), int((payload.get("usage") or {}).get("total_tokens") or 0)


def call_llm_json(prompt: str, config: dict[str, Any]) -> LLMResult:
    provider = config["llm"]["provider"]
    if provider == "disabled_mock_for_tests":
        payload = {
            "proposal_id": "mock_proposal_for_tests",
            "seed_id": "seed_001",
            "agent_role": "medicinal_chemist",
            "proposed_transformation_class": "small_substituent_scan",
            "edit_site_id": "seed_001_site_01",
            "transformation_description": "Add a methyl group at an unprotected peripheral atom.",
            "expected_effect": "probe hydrophobic occupancy",
            "expected_preserved_interactions": [],
            "expected_new_interactions": [],
            "medchem_risk_prediction": "low",
            "priority": 1,
            "requires_tool": True,
            "tool_name": "rdkit_transform",
            "tool_arguments_json": {},
            "rationale_summary": "Conservative test output for schema validation.",
            "confidence": 0.5,
            "reject_if_conditions_json": ["invalid_valence", "binding_mode_broken"],
        }
        return LLMResult(True, json.dumps(payload), "mock_for_tests_only", [])
    ok, status = provider_available(config)
    if not ok:
        return LLMResult(False, "", "provider_unavailable", [status])
    if provider not in {"llama_cpp_server", "vllm", "openai_compatible_local"}:
        return LLMResult(False, "", "provider_call_not_implemented", [provider])
    url = config["llm"]["base_url"].rstrip("/") + "/v1/chat/completions"
    warnings: list[str] = []
    last_status = ""
    for attempt in range(int(config["llm"].get("retry_count", 0)) + 1):
        try:
            response = requests.post(
                url,
                headers=_headers(config),
                json=_json_payload(prompt, config),
                timeout=int(config["llm"].get("timeout_seconds", 120)),
            )
        except requests.RequestException as exc:
            last_status = f"{provider}_chat_unreachable:{exc.__class__.__name__}"
            warnings.append(last_status)
            continue
        last_status = f"{provider}_chat_status_{response.status_code}"
        if not response.ok:
            warnings.append(f"{last_status}:{response.text[:300]}")
            continue
        raw_text, token_count = _extract_content(response.json())
        if raw_text:
            return LLMResult(True, raw_text, last_status, warnings, token_count)
        warnings.append(f"{last_status}:empty_content")
    return LLMResult(False, "", last_status or "llm_call_failed", warnings)
