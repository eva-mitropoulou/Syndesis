from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

from egfr_dockingforge.common.io import write_table


def score_md_stability(candidates: pd.DataFrame, metrics: pd.DataFrame, persistence: pd.DataFrame, paths: dict[str, Path]) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    cidx = candidates.set_index("md_candidate_id")
    pidx = persistence.set_index(["md_candidate_id", "md_system_id", "replicate_id"]) if not persistence.empty else {}
    rows = []
    for row in metrics.to_dict("records"):
        cand = cidx.loc[row["md_candidate_id"]]
        complete = row.get("trajectory_analysis_status") == "complete"
        if complete:
            ligand_component = max(0.0, 1.0 - float(row["ligand_rmsd_median_angstrom"]) / 6.0)
            pocket_component = max(0.0, 1.0 - float(row.get("protein_backbone_rmsd_median_angstrom") or 0.0) / 6.0)
            retention = float(row.get("fraction_frames_inside_pocket") or 0.0)
            key = pidx.loc[(row["md_candidate_id"], row["md_system_id"], row["replicate_id"])] if hasattr(pidx, "loc") and (row["md_candidate_id"], row["md_system_id"], row["replicate_id"]) in pidx.index else {}
            persistence_score = float(key.get("binding_mode_persistence_score", 0.0)) if hasattr(key, "get") else 0.0
            score = 0.35 * ligand_component + 0.25 * pocket_component + 0.25 * retention + 0.15 * persistence_score
            stable = (
                row["ligand_rmsd_median_angstrom"] <= 3.0
                and row["ligand_rmsd_p95_angstrom"] <= 5.0
                and retention >= 0.90
            )
            label = "md_stable" if stable else "md_unstable"
            reason = "" if stable else "trajectory_completed_but_stability_thresholds_not_met"
        else:
            ligand_component = 0.0
            pocket_component = 0.0
            persistence_score = 0.0
            retention = 0.0
            score = 0.0
            label = "md_failed_setup"
            reason = "missing_production_trajectory_or_analysis_failure"
        rows.append(
            {
                "md_candidate_id": row["md_candidate_id"],
                "molecule_id": cand["molecule_id"],
                "source": cand["source"],
                "best_pose_id": cand["best_pose_id"],
                "replicate_id": row["replicate_id"],
                "md_stability_label": label,
                "md_stability_score": score,
                "ligand_rmsd_component": ligand_component,
                "pocket_stability_component": pocket_component,
                "interaction_persistence_component": persistence_score,
                "pocket_retention_component": retention,
                "binding_mode_preserved_flag": label == "md_stable",
                "ligand_left_pocket_flag": row.get("ligand_left_pocket_flag"),
                "md_acceptance_flag": label == "md_stable",
                "md_rejection_reason": reason,
                "warnings_json": json.dumps([] if complete else ["no_md_stability_evidence_generated"]),
            }
    )
    labels = pd.DataFrame(rows)
    summaries = []
    has_ligand_rmsd = "ligand_rmsd_median_angstrom" in metrics.columns
    for cand in candidates.to_dict("records"):
        g = labels[labels["md_candidate_id"].eq(cand["md_candidate_id"])]
        cand_metrics = metrics[metrics["md_candidate_id"].eq(cand["md_candidate_id"])]
        cand_persistence = persistence[persistence["md_candidate_id"].eq(cand["md_candidate_id"])] if not persistence.empty else pd.DataFrame()
        median_ligand_rmsd = (
            float(cand_metrics["ligand_rmsd_median_angstrom"].median())
            if has_ligand_rmsd and not cand_metrics.empty
            else None
        )
        median_key_occupancy = (
            float(cand_persistence["key_interaction_occupancy_mean"].median())
            if "key_interaction_occupancy_mean" in cand_persistence.columns and not cand_persistence.empty
            else 0.0
        )
        median_binding_persistence = (
            float(cand_persistence["binding_mode_persistence_score"].median())
            if "binding_mode_persistence_score" in cand_persistence.columns and not cand_persistence.empty
            else float(g["interaction_persistence_component"].median()) if not g.empty else 0.0
        )
        summaries.append(
            {
                "md_candidate_id": cand["md_candidate_id"],
                "molecule_id": cand["molecule_id"],
                "source": cand["source"],
                "parent_molecule_id": cand["parent_molecule_id"],
                "num_replicates_completed": int(g["md_stability_label"].ne("md_failed_setup").sum()),
                "num_replicates_stable": int(g["md_stability_label"].eq("md_stable").sum()),
                "num_replicates_moderately_stable": 0,
                "num_replicates_unstable": int(g["md_stability_label"].eq("md_unstable").sum()),
                "replicate_success_fraction": float(g["md_stability_label"].ne("md_failed_setup").mean()) if not g.empty else 0.0,
                "median_md_stability_score": float(g["md_stability_score"].median()) if not g.empty else 0.0,
                "min_md_stability_score": float(g["md_stability_score"].min()) if not g.empty else 0.0,
                "median_ligand_rmsd": median_ligand_rmsd,
                "median_key_interaction_occupancy": median_key_occupancy,
                "median_binding_mode_persistence_score": median_binding_persistence,
                "final_md_label": g["md_stability_label"].iloc[0] if not g.empty else "md_failed_setup",
                "final_md_decision": "pass_md_stability" if not g.empty and g["md_stability_label"].iloc[0] == "md_stable" else "fail_md_stability",
                "final_md_reason": g["md_rejection_reason"].iloc[0] if not g.empty else "no_md_metrics",
                "recommended_for_final_dossier": bool(not g.empty and g["md_stability_label"].iloc[0] == "md_stable"),
                "warnings_json": json.dumps([] if not g.empty and g["md_stability_label"].iloc[0] != "md_failed_setup" else ["not_recommended_without_successful_md"]),
            }
        )
    summary = pd.DataFrame(summaries)
    summary_idx = summary.set_index("md_candidate_id") if not summary.empty else pd.DataFrame()
    post = pd.DataFrame(
        [
            {
                "analog_id": cand["molecule_id"],
                "md_candidate_id": cand["md_candidate_id"],
                "molecule_id": cand["molecule_id"],
                "parent_molecule_id": cand["parent_molecule_id"],
                "strategy_name": cand.get("stage10_strategy_name", ""),
                "source": cand["source"],
                "pre_md_accepted_flag": False,
                "final_md_label": summary_idx.loc[cand["md_candidate_id"], "final_md_label"] if cand["md_candidate_id"] in summary_idx.index else "md_failed_setup",
                "final_md_decision": summary_idx.loc[cand["md_candidate_id"], "final_md_decision"] if cand["md_candidate_id"] in summary_idx.index else "fail_md_stability",
                "md_stability_score": float(summary_idx.loc[cand["md_candidate_id"], "median_md_stability_score"]) if cand["md_candidate_id"] in summary_idx.index else 0.0,
                "binding_mode_persistence_score": float(summary_idx.loc[cand["md_candidate_id"], "median_binding_mode_persistence_score"]) if cand["md_candidate_id"] in summary_idx.index else 0.0,
                "key_interaction_occupancy_mean": float(summary_idx.loc[cand["md_candidate_id"], "median_key_interaction_occupancy"]) if cand["md_candidate_id"] in summary_idx.index else 0.0,
                "ligand_left_pocket_flag": None,
                "accepted_post_md_flag": bool(cand["md_candidate_id"] in summary_idx.index and summary_idx.loc[cand["md_candidate_id"], "final_md_decision"] == "pass_md_stability"),
                "post_md_rejection_reason": "" if cand["md_candidate_id"] in summary_idx.index and summary_idx.loc[cand["md_candidate_id"], "final_md_decision"] == "pass_md_stability" else "md_stability_thresholds_not_met_or_analysis_failed",
                "warnings_json": json.dumps([] if cand["md_candidate_id"] in summary_idx.index and summary_idx.loc[cand["md_candidate_id"], "final_md_decision"] == "pass_md_stability" else ["stage10B_input_records_md_not_accepted"]),
            }
            for cand in candidates.to_dict("records")
        ]
    )
    write_table(paths["processed"] / "md_pose_stability_labels.parquet", labels)
    write_table(paths["processed"] / "md_pose_stability_labels.csv", labels)
    write_table(paths["processed"] / "md_candidate_summary.parquet", summary)
    write_table(paths["processed"] / "md_candidate_summary.csv", summary)
    write_table(paths["processed"] / "stage10_post_md_input.parquet", post)
    write_table(paths["processed"] / "stage10_post_md_input.csv", post)
    return labels, summary, post
