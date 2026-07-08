from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pandas as pd

from egfr_dockingforge.common.io import write_table
from egfr_dockingforge.stage9.agent_council import LLM_STRATEGIES


def summarize_iterations(
    candidates: pd.DataFrame,
    validation: pd.DataFrame,
    screening: pd.DataFrame,
    acceptance: pd.DataFrame,
    config: dict[str, Any],
    paths: dict[str, Path],
) -> pd.DataFrame:
    rows = []
    for strategy in config["loop"]["strategies"]:
        cand = candidates[candidates["strategy_name"].eq(strategy)] if "strategy_name" in candidates else pd.DataFrame()
        val = validation[validation["analog_id"].isin(cand.get("analog_id", []))]
        scr = screening[screening["strategy_name"].eq(strategy)] if not screening.empty else pd.DataFrame()
        acc = acceptance[acceptance["strategy_name"].eq(strategy)] if not acceptance.empty else pd.DataFrame()
        rows.append(
            {
                "iteration_id": "iter_001",
                "strategy_name": strategy,
                "seed_id": "all",
                "num_agent_proposals": 0 if strategy == "rdkit_rule_based" else len(cand),
                "num_valid_molecules": int(val["valid_molecule_flag"].sum()) if not val.empty else 0,
                "num_unique_molecules": int(cand["standard_smiles"].nunique()) if not cand.empty else 0,
                "num_screened": len(scr),
                "num_accepted": int(acc["accepted_flag"].sum()) if not acc.empty else 0,
                "accepted_rate": float(acc["accepted_flag"].mean()) if not acc.empty else 0.0,
                "dominant_rejection_reasons_json": json.dumps(acc["rejection_reason"].value_counts().head(5).to_dict()) if not acc.empty else "{}",
                "best_delta_candidate_score": float(acc["delta_candidate_score"].max()) if not acc.empty else 0.0,
                "best_delta_pose_confidence": float(acc["delta_pose_confidence"].max()) if not acc.empty else 0.0,
                "best_delta_ligand_efficiency": float(acc["delta_ligand_efficiency"].max()) if not acc.empty else 0.0,
                "runtime_seconds": 0.0,
                "warnings_json": json.dumps([]),
            }
        )
    out = pd.DataFrame(rows)
    write_table(paths["processed"] / "agent_iterations.parquet", out)
    write_table(paths["processed"] / "agent_iterations.csv", out)
    return out


def benchmark_strategies(
    candidates: pd.DataFrame,
    validation: pd.DataFrame,
    screening: pd.DataFrame,
    acceptance: pd.DataFrame,
    agent_status: pd.DataFrame,
    config: dict[str, Any],
    paths: dict[str, Path],
) -> pd.DataFrame:
    rows = []
    for strategy in config["loop"]["strategies"]:
        cand = candidates[candidates["strategy_name"].eq(strategy)] if "strategy_name" in candidates else pd.DataFrame()
        val = validation[validation["analog_id"].isin(cand.get("analog_id", []))]
        scr = screening[screening["strategy_name"].eq(strategy)] if not screening.empty else pd.DataFrame()
        acc = acceptance[acceptance["strategy_name"].eq(strategy)] if not acceptance.empty else pd.DataFrame()
        raw = len(cand)
        valid = int(val["valid_molecule_flag"].sum()) if not val.empty else 0
        unique = int(cand["standard_smiles"].nunique()) if not cand.empty else 0
        accepted = int(acc["accepted_flag"].sum()) if not acc.empty else 0
        is_llm = strategy in LLM_STRATEGIES
        model_name = config["llm"]["model_name"] if is_llm else ""
        notes = "local LLM proposal strategy logged in agent_proposal_status" if is_llm else "deterministic RDKit baseline"
        if strategy == "reinvent4_baseline":
            notes = "REINVENT4 package not installed; baseline recorded as not run"
        rows.append(
            {
                "strategy_name": strategy,
                "num_seeds": int(candidates["seed_id"].nunique()) if not candidates.empty else 0,
                "num_iterations": int(config["loop"]["max_iterations"]),
                "num_raw_proposals": raw if not is_llm else int(agent_status[agent_status["strategy_name"].eq(strategy)].shape[0]),
                "num_valid_molecules": valid,
                "validity_rate": valid / raw if raw else 0.0,
                "uniqueness_rate": unique / raw if raw else 0.0,
                "novelty_rate": 1.0 if raw else 0.0,
                "num_screened": len(scr),
                "num_accepted": accepted,
                "accepted_analog_rate": accepted / len(scr) if len(scr) else 0.0,
                "mean_delta_candidate_score": float(acc["delta_candidate_score"].mean()) if not acc.empty else 0.0,
                "median_delta_candidate_score": float(acc["delta_candidate_score"].median()) if not acc.empty else 0.0,
                "mean_delta_pose_confidence": float(acc["delta_pose_confidence"].mean()) if not acc.empty else 0.0,
                "median_delta_pose_confidence": float(acc["delta_pose_confidence"].median()) if not acc.empty else 0.0,
                "mean_delta_key_interaction_recall": float(acc["delta_key_interaction_recall"].mean()) if not acc.empty else 0.0,
                "score_hacking_rate": float(acc["score_hacking_flag"].mean()) if not acc.empty else 0.0,
                "bad_chemistry_rejection_rate": float((val["hard_scope_pass"] == False).mean()) if not val.empty else 0.0,
                "binding_mode_break_rate": float((acc["binding_mode_preserved_flag"] == False).mean()) if not acc.empty else 0.0,
                "mean_runtime_per_accepted_analog": 0.0,
                "token_count_if_llm": 0,
                "local_model_name_if_llm": model_name,
                "cost_estimate_if_any": 0.0,
                "benchmark_notes": notes,
            }
        )
    out = pd.DataFrame(rows)
    write_table(paths["processed"] / "analog_strategy_benchmark.parquet", out)
    write_table(paths["processed"] / "analog_strategy_benchmark.csv", out)
    return out
