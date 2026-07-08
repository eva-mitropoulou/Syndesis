from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pandas as pd

from egfr_dockingforge.common.io import write_table
from egfr_dockingforge.stage10.strategy_normalization import STRATEGIES


def build_ablation_manifest(inputs: dict[str, pd.DataFrame | None], config: dict[str, Any], paths: dict[str, Path]) -> pd.DataFrame:
    seeds = inputs["stage9_seeds"]
    md_available = inputs.get("stage11_md_metrics") is not None
    rows = []
    for i, name in enumerate(STRATEGIES, start=1):
        uses_llm = name in {"single_agent", "council_loop", "council_plus_prolif", "council_plus_prolif_pose_confidence", "council_plus_prolif_pose_confidence_md"}
        uses_council = name.startswith("council")
        md = name.endswith("_md")
        rows.append(
            {
                "strategy_id": f"strat_{i:02d}",
                "strategy_name": name,
                "strategy_family": "agentic" if uses_llm else "baseline",
                "uses_llm": uses_llm,
                "uses_council": uses_council,
                "uses_tool_feedback": name in {"council_plus_prolif", "council_plus_prolif_pose_confidence", "council_plus_prolif_pose_confidence_md"},
                "uses_rdkit_validation": True,
                "uses_docking_score": name != "random_analog_enumeration",
                "uses_gnina": name not in {"random_analog_enumeration", "docking_score_only_optimization"},
                "uses_prolif_constraint": name in {"council_plus_prolif", "council_plus_prolif_pose_confidence", "council_plus_prolif_pose_confidence_md", "rdkit_rule_based"},
                "uses_pose_confidence": name in {"council_plus_prolif_pose_confidence", "council_plus_prolif_pose_confidence_md", "rdkit_rule_based"},
                "uses_md_filter": md,
                "uses_medchem_filter": True,
                "generation_budget": int(config["benchmark"]["generation_budget"]),
                "screening_budget": int(config["benchmark"]["screening_budget"]),
                "llm_model_name": config.get("llm_model_name", ""),
                "seed_set_id": config["benchmark"]["seed_set_id"],
                "num_seed_scaffolds": len(seeds) if seeds is not None else 0,
                "num_iterations": 1,
                "enabled_flag": False if md and not md_available else True,
                "md_status": "available" if md_available else ("pending_stage11" if md else "not_applicable"),
                "notes": "pre-MD benchmark" if not md else "post-MD ablation pending Stage 11",
            }
        )
    out = pd.DataFrame(rows)
    write_table(paths["processed"] / "ablation_strategy_manifest.parquet", out)
    write_table(paths["processed"] / "ablation_strategy_manifest.csv", out)
    return out


def build_budget_audit(master: pd.DataFrame, manifest: pd.DataFrame, inputs: dict[str, pd.DataFrame | None], config: dict[str, Any], paths: dict[str, Path]) -> pd.DataFrame:
    seeds = inputs["stage9_seeds"]
    rows = []
    for strat in manifest.to_dict("records"):
        for seed_id in seeds["seed_id"].tolist():
            g = master[(master["strategy_id"] == strat["strategy_id"]) & (master["seed_id"] == seed_id)]
            valid_unique = g[g["valid_molecule_flag"] & g["unique_flag"]]
            screened = g[g["best_pose_confidence"].notna()]
            violation = ""
            if len(g) < strat["generation_budget"] and strat["enabled_flag"]:
                violation = "generated_fewer_than_budget"
            rows.append(
                {
                    "strategy_id": strat["strategy_id"],
                    "seed_id": seed_id,
                    "num_raw_proposals": len(g),
                    "num_valid_unique_analogs": len(valid_unique),
                    "num_screened_analogs": len(screened),
                    "num_docking_tasks": len(screened),
                    "num_gnina_tasks": len(screened),
                    "num_prolif_tasks": len(screened),
                    "num_pose_confidence_predictions": len(screened),
                    "num_md_tasks_if_available": 0,
                    "llm_token_count": 0,
                    "walltime_seconds": 0.0,
                    "gpu_seconds": 0.0,
                    "cpu_seconds": 0.0,
                    "budget_normalized_flag": len(g) <= strat["generation_budget"] and len(screened) <= strat["screening_budget"],
                    "budget_violation_reason": violation,
                }
            )
    out = pd.DataFrame(rows)
    write_table(paths["processed"] / "strategy_budget_audit.parquet", out)
    write_table(paths["processed"] / "strategy_budget_audit.csv", out)
    return out
