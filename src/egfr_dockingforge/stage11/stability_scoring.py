from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

from egfr_dockingforge.common.io import write_table


def score_md_stability(candidates: pd.DataFrame, metrics: pd.DataFrame, persistence: pd.DataFrame, paths: dict[str, Path], config: dict | None = None) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    cidx = candidates.set_index("md_candidate_id")
    pidx = persistence.set_index(["md_candidate_id", "md_system_id", "replicate_id"]) if not persistence.empty else {}
    # Stability thresholds (config-driven). The interaction thresholds default to
    # 0.0 when no config is passed, which preserves the old RMSD-only behaviour for
    # any legacy caller; the Stage-11 CLI passes config so the interaction gate is active.
    stab = (config or {}).get("stability", {}) if config else {}
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
            hinge_occ = float(key.get("hinge_interaction_occupancy", 0.0)) if hasattr(key, "get") else 0.0
            key_occ = float(key.get("key_interaction_occupancy_mean", 0.0)) if hasattr(key, "get") else 0.0
            score = 0.35 * ligand_component + 0.25 * pocket_component + 0.25 * retention + 0.15 * persistence_score
            # Stability now requires BOTH pose geometry (RMSD + pocket retention) AND
            # maintenance of the defining EGFR interactions. Without the interaction
            # requirement a pose could sit low-RMSD in the pocket while having lost the
            # hinge H-bond entirely (observed: a control passed at hinge occupancy
            # ~0.01) -- a geometrically-stable but mechanistically-wrong "binder".
            # Thresholds are config-driven (stability.*), defaulting to lenient values
            # that only fail a genuinely lost hinge/key-contact pose.
            rmsd_ok = (
                row["ligand_rmsd_median_angstrom"] <= float(stab.get("ligand_rmsd_median_stable_angstrom", 3.0))
                and row["ligand_rmsd_p95_angstrom"] <= float(stab.get("ligand_rmsd_p95_stable_angstrom", 5.0))
                and retention >= float(stab.get("fraction_frames_inside_pocket_stable", 0.90))
            )
            # Interaction thresholds default to 0.0 (gate disabled) when no config is
            # supplied, so legacy RMSD-only behaviour is preserved; the Stage-11 config
            # sets 0.30 / 0.50 to require a maintained hinge + key-contact set.
            hinge_thr = float(stab.get("hinge_interaction_occupancy_stable", 0.0))
            key_thr = float(stab.get("key_interaction_occupancy_mean_stable", 0.0))
            interaction_ok = (hinge_occ >= hinge_thr and key_occ >= key_thr)
            stable = bool(rmsd_ok and interaction_ok)
            label = "md_stable" if stable else "md_unstable"
            if stable:
                reason = ""
            elif not rmsd_ok:
                reason = "pose_geometry_unstable_rmsd_or_pocket"
            else:
                reason = f"binding_interactions_not_maintained_hinge{hinge_occ:.2f}_key{key_occ:.2f}"
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
    # An empty cohort (no MD metrics reached scoring — e.g. every system failed
    # equilibration) must NOT KeyError on labels["md_candidate_id"]; guarantee the
    # column exists so downstream .eq() filters return empty rather than crashing.
    if labels.empty:
        labels = pd.DataFrame(columns=["md_candidate_id", "molecule_id", "md_stability_label", "interaction_persistence_component"])
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
        n_completed = int(g["md_stability_label"].ne("md_failed_setup").sum())
        n_stable = int(g["md_stability_label"].eq("md_stable").sum())
        n_unstable = int(g["md_stability_label"].eq("md_unstable").sum())
        # Final verdict = MAJORITY of completed replicates stable. A single-replicate
        # verdict (the old `iloc[0]`) was order-dependent and arbitrary: it labelled
        # a 2/3-stable finalist by whichever replicate happened to sort first. The
        # majority rule is the standard robustness criterion for independent
        # replicates and is calibrated by the controls (all 3 known EGFR inhibitors
        # pass it). Ties / all-failed -> not stable.
        if n_completed == 0:
            final_label = "md_failed_setup"
            final_reason = "no_completed_replicates"
        elif n_stable * 2 > n_completed:
            final_label = "md_stable"
            final_reason = ""
        else:
            final_label = "md_unstable"
            final_reason = f"stable_in_{n_stable}_of_{n_completed}_replicates_below_majority"
        final_decision = "pass_md_stability" if final_label == "md_stable" else "fail_md_stability"
        summaries.append(
            {
                "md_candidate_id": cand["md_candidate_id"],
                "molecule_id": cand["molecule_id"],
                "source": cand["source"],
                "parent_molecule_id": cand["parent_molecule_id"],
                "num_replicates_completed": n_completed,
                "num_replicates_stable": n_stable,
                "num_replicates_moderately_stable": 0,
                "num_replicates_unstable": n_unstable,
                "replicate_success_fraction": float(g["md_stability_label"].ne("md_failed_setup").mean()) if not g.empty else 0.0,
                "median_md_stability_score": float(g["md_stability_score"].median()) if not g.empty else 0.0,
                "min_md_stability_score": float(g["md_stability_score"].min()) if not g.empty else 0.0,
                "median_ligand_rmsd": median_ligand_rmsd,
                "median_key_interaction_occupancy": median_key_occupancy,
                "median_binding_mode_persistence_score": median_binding_persistence,
                "final_md_label": final_label,
                "final_md_decision": final_decision,
                "final_md_reason": final_reason,
                "recommended_for_final_dossier": bool(final_label == "md_stable"),
                "warnings_json": json.dumps([] if final_label != "md_failed_setup" else ["not_recommended_without_successful_md"]),
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
