from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pandas as pd

from egfr_dockingforge.common.io import ensure_dir, write_table
from egfr_dockingforge.stage9.agent_prompts import AGENT_ROLES, compact_prompt
from egfr_dockingforge.stage9.agent_schemas import parse_agent_json
from egfr_dockingforge.stage9.llm_clients import call_llm_json, provider_available


LLM_STRATEGIES = {
    "single_agent": ["medicinal_chemist"],
    "council_no_feedback": ["medicinal_chemist", "docking_scientist", "interaction_analyst", "admet_safety_critic", "skeptical_reviewer"],
    "council_with_tool_feedback": ["medicinal_chemist", "docking_scientist", "interaction_analyst", "admet_safety_critic", "skeptical_reviewer"],
    "council_interaction_constrained": ["medicinal_chemist", "docking_scientist", "interaction_analyst", "admet_safety_critic", "skeptical_reviewer", "orchestrator"],
}


def run_agent_loop(seeds: pd.DataFrame, edit_sites: pd.DataFrame, config: dict[str, Any], paths: dict[str, Path]) -> pd.DataFrame:
    records = []
    out_jsonl = paths["processed"] / "analog_proposals.jsonl"
    ensure_dir(out_jsonl.parent)
    ok, status = provider_available(config)
    configured_llm_strategies = [strategy for strategy in config["loop"]["strategies"] if strategy in LLM_STRATEGIES]
    if configured_llm_strategies and not ok:
        raise RuntimeError(f"Stage 9 local LLM provider is unavailable for configured strategies {configured_llm_strategies}: {status}")
    with out_jsonl.open("w", encoding="utf-8") as handle:
        for strategy, roles in LLM_STRATEGIES.items():
            if strategy not in config["loop"]["strategies"]:
                continue
            for seed in seeds.to_dict("records"):
                site_rows = edit_sites[edit_sites["seed_id"].eq(seed["seed_id"])].head(3).to_dict("records")
                for role in roles:
                    prompt = compact_prompt(role, seed, site_rows)
                    result = call_llm_json(prompt, config)
                    parsed = parse_agent_json(result.raw_text)
                    record = {
                        "proposal_id": f"{strategy}_{seed['seed_id']}_{role}",
                        "iteration_id": "iter_001",
                        "strategy_name": strategy,
                        "seed_id": seed["seed_id"],
                        "agent_role": role,
                        "raw_agent_output": result.raw_text,
                        "parsed_json": json.dumps(parsed.parsed_json or {}),
                        "schema_valid": parsed.schema_valid,
                        "repair_attempts": parsed.repair_attempts,
                        "proposal_status": parsed.proposal_status if result.ok else result.status,
                        "warnings_json": json.dumps(result.warnings or [status] if not ok else result.warnings),
                        "token_count": result.token_count,
                    }
                    handle.write(json.dumps(record, sort_keys=True) + "\n")
                    records.append(record)
    out = pd.DataFrame(records)
    write_table(paths["processed"] / "agent_proposal_status.parquet", out)
    write_table(paths["processed"] / "agent_proposal_status.csv", out)
    roles = pd.DataFrame([{"agent_role": k, "definition": v} for k, v in AGENT_ROLES.items()])
    write_table(paths["processed"] / "agent_roles.parquet", roles)
    return out
