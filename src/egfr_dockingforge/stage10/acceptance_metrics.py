from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pandas as pd

from egfr_dockingforge.common.io import write_table
from egfr_dockingforge.stage10.metric_definitions import accepted_pre_md, score_hacking
from egfr_dockingforge.stage10.strategy_normalization import normalize_strategy_name


def build_master_table(inputs: dict[str, pd.DataFrame | None], manifest: pd.DataFrame, config: dict[str, Any], paths: dict[str, Path]) -> pd.DataFrame:
    candidates = inputs["stage9_candidates"].copy()
    validation = inputs["stage9_validation"].copy()
    screening = inputs["stage9_screening"].copy()
    acceptance = inputs["stage9_acceptance"].copy()
    seeds = inputs["stage9_seeds"].copy()
    df = candidates.merge(validation, on="analog_id", how="left").merge(screening, on=["analog_id","seed_id","strategy_name","iteration_id"], how="left").merge(acceptance, on=["analog_id","seed_id","strategy_name","iteration_id"], how="left", suffixes=("", "_acc"))
    seed_idx = seeds.set_index("seed_id")
    strat_id = manifest.set_index("strategy_name")["strategy_id"].to_dict()
    rows = []
    for row in df.to_dict("records"):
        norm = normalize_strategy_name(row["strategy_name"])
        seed = seed_idx.loc[row["seed_id"]]
        item = {
            "analog_id": row["analog_id"],
            "strategy_id": strat_id.get(norm, ""),
            "strategy_name": norm,
            "seed_id": row["seed_id"],
            "parent_molecule_id": row["parent_molecule_id"],
            "parent_smiles": row["parent_smiles"],
            "analog_smiles": row["analog_smiles"],
            "standard_smiles": row["standard_smiles"],
            "inchi_key": row["inchi_key"],
            "scaffold_id": seed.get("scaffold_id", row["seed_id"]),
            "valid_molecule_flag": bool(row.get("valid_molecule_flag", False)),
            "unique_flag": row.get("uniqueness_status", "") in {"unique", ""},
            "hard_scope_pass": bool(row.get("hard_scope_pass", False)),
            "covalent_warhead_flag": bool(row.get("covalent_warhead_flag", False)),
            "reactive_flag": bool(row.get("reactive_flag", False)),
            "pains_flag": bool(row.get("pains_flag", False)),
            "property_window_pass": bool(row.get("property_window_pass", False)),
            "medchem_risk_score": float(row.get("medchem_risk_score", 1.0)),
            "parent_tanimoto": float(row.get("parent_tanimoto", 0.0)),
            "closest_known_egfr_ligand": row.get("closest_known_egfr_ligand", ""),
            "novelty_bucket": row.get("novelty_status", "stage9_analog"),
            "best_docking_score": row.get("best_docking_score"),
            "best_gnina_cnnscore": row.get("best_gnina_cnnscore"),
            "best_gnina_cnnaffinity": row.get("best_gnina_cnnaffinity"),
            "best_pose_confidence": row.get("best_pose_confidence"),
            "best_key_interaction_recall_consensus": row.get("best_key_interaction_recall_consensus"),
            "best_ifp_tanimoto_to_consensus": row.get("best_ifp_tanimoto_to_consensus"),
            "binding_mode_preserved_flag": bool(row.get("binding_mode_preserved_flag", False)),
            "ligand_efficiency": row.get("ligand_efficiency"),
            "parent_candidate_score": row.get("parent_candidate_score"),
            "analog_candidate_score": row.get("analog_candidate_score"),
            "delta_candidate_score": row.get("delta_candidate_score", 0.0),
            "delta_pose_confidence": row.get("delta_pose_confidence", 0.0),
            "delta_gnina_cnnscore": row.get("delta_gnina_cnnscore", 0.0),
            "delta_key_interaction_recall": row.get("delta_key_interaction_recall", 0.0),
            "delta_ligand_efficiency": row.get("delta_ligand_efficiency", 0.0),
            "md_stability_label_if_available": None,
            "md_key_interaction_persistence_if_available": None,
            "accepted_post_md_flag": None,
            "rejection_reason": row.get("rejection_reason", ""),
            "warnings_json": json.dumps([]),
        }
        item["accepted_pre_md_flag"] = accepted_pre_md(
            item,
            float(config["acceptance"]["min_pose_confidence"]),
            float(config["acceptance"]["max_medchem_risk"]),
            float(config["acceptance"].get("min_delta_candidate_score", -0.02)),
            float(config["acceptance"].get("min_delta_ligand_efficiency", -0.10)),
        )
        item["score_hacking_flag"] = score_hacking(item)
        rows.append(item)
    out = pd.DataFrame(rows)
    write_table(paths["processed"] / "analog_benchmark_master.parquet", out)
    write_table(paths["processed"] / "analog_benchmark_master.csv", out)
    return out


def compute_seed_and_strategy_metrics(master: pd.DataFrame, manifest: pd.DataFrame, paths: dict[str, Path]) -> tuple[pd.DataFrame, pd.DataFrame]:
    seed_rows = []
    seeds = sorted(master["seed_id"].dropna().unique().tolist())
    for strat in manifest.to_dict("records"):
        for seed_id in seeds:
            g = master[(master["seed_id"].eq(seed_id)) & (master["strategy_id"].eq(strat["strategy_id"]))]
            accepted = g[g["accepted_pre_md_flag"]] if not g.empty else pd.DataFrame()
            seed_rows.append({
                "seed_id": seed_id, "strategy_id": strat["strategy_id"], "strategy_name": strat["strategy_name"],
                "num_raw_proposals": len(g), "num_valid_molecules": int(g["valid_molecule_flag"].sum()) if not g.empty else 0,
                "num_unique_molecules": int(g["unique_flag"].sum()) if not g.empty else 0, "num_screened": int(g["best_pose_confidence"].notna().sum()) if not g.empty else 0,
                "num_accepted_pre_md": int(g["accepted_pre_md_flag"].sum()) if not g.empty else 0, "num_accepted_post_md_if_available": None,
                "accepted_rate_pre_md": float(g["accepted_pre_md_flag"].sum()/max(int(g["valid_molecule_flag"].sum()),1)) if not g.empty else 0.0,
                "accepted_rate_post_md_if_available": None, "score_hacking_rate": float(g["score_hacking_flag"].mean()) if not g.empty else 0.0,
                "best_delta_candidate_score": float(g["delta_candidate_score"].max()) if not g.empty else 0.0,
                "best_delta_pose_confidence": float(g["delta_pose_confidence"].max()) if not g.empty else 0.0,
                "best_delta_key_interaction_recall": float(g["delta_key_interaction_recall"].max()) if not g.empty else 0.0,
                "best_delta_ligand_efficiency": float(g["delta_ligand_efficiency"].max()) if not g.empty else 0.0,
                "best_accepted_analog_id": accepted["analog_id"].iloc[0] if not accepted.empty else "",
                "dominant_rejection_reason": g["rejection_reason"].fillna("").value_counts().index[0] if not g.empty else "no_analogs_generated",
                "warnings_json": json.dumps([] if not g.empty else ["strategy_generated_no_analogs"]),
            })
    seed_metrics = pd.DataFrame(seed_rows)
    strategy_rows = []
    for strat in manifest.to_dict("records"):
        g = master[master["strategy_id"].eq(strat["strategy_id"])]
        valid = int(g["valid_molecule_flag"].sum()) if not g.empty else 0
        unique = int(g["unique_flag"].sum()) if not g.empty else 0
        accepted = int(g["accepted_pre_md_flag"].sum()) if not g.empty else 0
        strategy_rows.append({
            "strategy_id": strat["strategy_id"], "strategy_name": strat["strategy_name"],
            "num_seeds": strat["num_seed_scaffolds"], "num_raw_proposals": len(g), "num_valid_molecules": valid,
            "num_unique_molecules": unique, "num_screened": int(g["best_pose_confidence"].notna().sum()) if not g.empty else 0,
            "num_pre_md_accepted": accepted, "num_post_md_accepted_if_available": None,
            "validity_rate": valid/max(len(g),1), "uniqueness_rate": unique/max(len(g),1), "novelty_rate": 1.0 if len(g) else 0.0,
            "diversity_score": float(1 - g["parent_tanimoto"].mean()) if len(g) else 0.0,
            "accepted_analog_rate_pre_md": accepted/max(unique,1),
            "accepted_analog_rate_post_md_if_available": None,
            "score_hacking_rate": float(g["score_hacking_flag"].mean()) if len(g) else 0.0,
            "binding_mode_break_rate": float((g["binding_mode_preserved_flag"] == False).mean()) if len(g) else 0.0,
            "bad_chemistry_rejection_rate": float((g["hard_scope_pass"] == False).mean()) if len(g) else 0.0,
            "medchem_risk_rejection_rate": float((g["medchem_risk_score"] > 0.4).mean()) if len(g) else 0.0,
            "mean_delta_candidate_score": float(g["delta_candidate_score"].mean()) if len(g) else 0.0,
            "median_delta_candidate_score": float(g["delta_candidate_score"].median()) if len(g) else 0.0,
            "mean_delta_pose_confidence": float(g["delta_pose_confidence"].mean()) if len(g) else 0.0,
            "median_delta_pose_confidence": float(g["delta_pose_confidence"].median()) if len(g) else 0.0,
            "mean_delta_key_interaction_recall": float(g["delta_key_interaction_recall"].mean()) if len(g) else 0.0,
            "median_delta_key_interaction_recall": float(g["delta_key_interaction_recall"].median()) if len(g) else 0.0,
            "mean_delta_ligand_efficiency": float(g["delta_ligand_efficiency"].mean()) if len(g) else 0.0,
            "median_delta_ligand_efficiency": float(g["delta_ligand_efficiency"].median()) if len(g) else 0.0,
            "mean_runtime_per_valid_analog": 0.0, "mean_runtime_per_accepted_analog": 0.0,
            "gpu_seconds_per_accepted_analog": 0.0, "llm_tokens_per_accepted_analog": 0.0,
            "md_status": strat["md_status"], "notes": strat["notes"],
        })
    strategy_metrics = pd.DataFrame(strategy_rows)
    write_table(paths["processed"] / "seed_strategy_metrics.parquet", seed_metrics)
    write_table(paths["processed"] / "seed_strategy_metrics.csv", seed_metrics)
    write_table(paths["processed"] / "strategy_metrics.parquet", strategy_metrics)
    write_table(paths["processed"] / "strategy_metrics.csv", strategy_metrics)
    return seed_metrics, strategy_metrics
