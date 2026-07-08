from __future__ import annotations

from pathlib import Path
from typing import Any

import pandas as pd

from egfr_dockingforge.common.io import write_json
from egfr_dockingforge.stage9.agent_council import run_agent_loop
from egfr_dockingforge.stage9.analog_acceptance import score_analog_acceptance
from egfr_dockingforge.stage9.analog_deduplication import deduplicate_analogs
from egfr_dockingforge.stage9.analog_tool_dispatch import dispatch_agent_proposals_to_molecules
from egfr_dockingforge.stage9.analog_screening_bridge import screen_analog_batch
from egfr_dockingforge.stage9.analog_validation import validate_analog_batch
from egfr_dockingforge.stage9.baseline_benchmark import benchmark_strategies, summarize_iterations
from egfr_dockingforge.stage9.edit_site_detection import detect_edit_sites
from egfr_dockingforge.stage9.load_stage_inputs import load_stage9_config, load_stage9_inputs, stage9_paths
from egfr_dockingforge.stage9.rdkit_transformations import enumerate_rule_based_analogs, write_transformation_library
from egfr_dockingforge.stage9.reinvent4_runner import run_reinvent4_baseline
from egfr_dockingforge.stage9.report_stage9 import write_stage9_report
from egfr_dockingforge.stage9.seed_selection import select_seed_scaffolds


def _load(config_path: str | Path) -> tuple[dict[str, Any], dict[str, Path], dict[str, pd.DataFrame]]:
    config = load_stage9_config(config_path)
    paths = stage9_paths(config)
    inputs = load_stage9_inputs(config)
    return config, paths, inputs


def select_analog_seeds_cli(config_path: str | Path) -> dict:
    config, paths, inputs = _load(config_path)
    seeds = select_seed_scaffolds(inputs, config, paths)
    return {"status": "complete", "seeds": len(seeds)}


def detect_edit_sites_cli(config_path: str | Path) -> dict:
    config, paths, _ = _load(config_path)
    seeds = pd.read_parquet(paths["processed"] / "analog_seed_scaffolds.parquet")
    edits = detect_edit_sites(seeds, config, paths)
    return {"status": "complete", "edit_sites": len(edits)}


def enumerate_rule_based_analogs_cli(config_path: str | Path) -> dict:
    config, paths, _ = _load(config_path)
    seeds = pd.read_parquet(paths["processed"] / "analog_seed_scaffolds.parquet")
    edits = pd.read_parquet(paths["processed"] / "edit_sites.parquet")
    write_transformation_library(paths)
    analogs = enumerate_rule_based_analogs(seeds, edits, config, paths)
    analogs = deduplicate_analogs(analogs)
    analogs.to_parquet(paths["processed"] / "analog_candidates.parquet", index=False)
    analogs.to_csv(paths["processed"] / "analog_candidates.csv", index=False)
    return {"status": "complete", "analogs": len(analogs)}


def run_reinvent4_baseline_cli(config_path: str | Path) -> dict:
    config, paths, _ = _load(config_path)
    status = run_reinvent4_baseline(config, paths)
    return {"status": "complete", "rows": len(status)}


def run_agent_analog_loop_cli(config_path: str | Path) -> dict:
    config, paths, _ = _load(config_path)
    seeds = pd.read_parquet(paths["processed"] / "analog_seed_scaffolds.parquet")
    edits = pd.read_parquet(paths["processed"] / "edit_sites.parquet")
    status = run_agent_loop(seeds, edits, config, paths)
    # Bridge: turn schema-valid agent proposals into concrete molecules and merge
    # them into the candidate pool so the agentic strategies are screened head to
    # head with the rule-based baseline (instead of producing zero molecules).
    agent_analogs = dispatch_agent_proposals_to_molecules(status, seeds, edits, config, paths)
    merged = _merge_into_candidate_pool(agent_analogs, paths)
    return {
        "status": "complete",
        "agent_rows": len(status),
        "agent_generated_analogs": int(len(agent_analogs)),
        "total_candidate_analogs": int(len(merged)),
    }


def _merge_into_candidate_pool(new_analogs: pd.DataFrame, paths: dict[str, Path]) -> pd.DataFrame:
    """Append newly-generated analogs to analog_candidates.parquet, de-duplicating
    on (seed_id, standard_smiles, strategy_name) so each strategy keeps its own
    molecules for a fair per-strategy accepted-rate comparison."""
    pool_path = paths["processed"] / "analog_candidates.parquet"
    existing = pd.read_parquet(pool_path) if pool_path.exists() else pd.DataFrame()
    if new_analogs.empty:
        return existing
    combined = pd.concat([existing, new_analogs], ignore_index=True) if not existing.empty else new_analogs
    combined = combined.drop_duplicates(["seed_id", "standard_smiles", "strategy_name"]).reset_index(drop=True)
    combined.to_parquet(pool_path, index=False)
    combined.to_csv(paths["processed"] / "analog_candidates.csv", index=False)
    return combined


def validate_analog_batch_cli(config_path: str | Path) -> dict:
    config, paths, _ = _load(config_path)
    candidates = pd.read_parquet(paths["processed"] / "analog_candidates.parquet")
    validation = validate_analog_batch(candidates, config, paths)
    return {"status": "complete", "validation_rows": len(validation)}


def screen_analog_batch_cli(config_path: str | Path) -> dict:
    config, paths, _ = _load(config_path)
    candidates = pd.read_parquet(paths["processed"] / "analog_candidates.parquet")
    validation = pd.read_parquet(paths["processed"] / "analog_validation.parquet")
    seeds = pd.read_parquet(paths["processed"] / "analog_seed_scaffolds.parquet")
    screening = screen_analog_batch(candidates, validation, seeds, config, paths)
    return {"status": "complete", "screened": len(screening)}


def score_analog_acceptance_cli(config_path: str | Path) -> dict:
    config, paths, _ = _load(config_path)
    screening = pd.read_parquet(paths["processed"] / "analog_screening_results.parquet")
    validation = pd.read_parquet(paths["processed"] / "analog_validation.parquet")
    seeds = pd.read_parquet(paths["processed"] / "analog_seed_scaffolds.parquet")
    acceptance = score_analog_acceptance(screening, validation, seeds, config, paths)
    return {"status": "complete", "acceptance_rows": len(acceptance), "accepted": int(acceptance["accepted_flag"].sum())}


def benchmark_analog_strategies_cli(config_path: str | Path) -> dict:
    config, paths, _ = _load(config_path)
    candidates = pd.read_parquet(paths["processed"] / "analog_candidates.parquet")
    validation = pd.read_parquet(paths["processed"] / "analog_validation.parquet")
    screening = pd.read_parquet(paths["processed"] / "analog_screening_results.parquet")
    acceptance = pd.read_parquet(paths["processed"] / "analog_acceptance.parquet")
    agent_status = pd.read_parquet(paths["processed"] / "agent_proposal_status.parquet")
    summarize_iterations(candidates, validation, screening, acceptance, config, paths)
    bench = benchmark_strategies(candidates, validation, screening, acceptance, agent_status, config, paths)
    return {"status": "complete", "strategies": len(bench)}


def report_stage9_cli(config_path: str | Path) -> dict:
    config, paths, _ = _load(config_path)
    report = write_stage9_report(paths)
    return {"status": "complete", "report": str(report)}


def run_stage9_all(config_path: str | Path) -> dict:
    select_analog_seeds_cli(config_path)
    detect_edit_sites_cli(config_path)
    enumerate_rule_based_analogs_cli(config_path)
    run_reinvent4_baseline_cli(config_path)
    run_agent_analog_loop_cli(config_path)
    validate_analog_batch_cli(config_path)
    screen_analog_batch_cli(config_path)
    score_analog_acceptance_cli(config_path)
    benchmark_analog_strategies_cli(config_path)
    summary = report_stage9_cli(config_path)
    config, paths, _ = _load(config_path)
    write_json(paths["processed"] / "stage9_summary.json", summary)
    return summary
