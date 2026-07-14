from __future__ import annotations

import json
import pandas as pd


def pose_decisions(conf: pd.DataFrame, normalized: pd.DataFrame, manifest: pd.DataFrame, config: dict, paths: dict) -> pd.DataFrame:
    df = conf.merge(normalized, on=["screening_pose_id", "molecule_id", "target_receptor_id", "receptor_state"], how="left")
    labels = []
    for row in df.to_dict("records"):
        if row["sanity_status"] != "pass":
            label = "physically_invalid"
        elif row["pose_confidence_probability"] >= config["triage"]["confidence_threshold"] and row["key_interaction_recall_consensus"] >= config["triage"]["interaction_f1_threshold"]:
            label = "accepted_pose_consistent"
        elif row["pose_confidence_probability"] < config["triage"]["confidence_threshold"]:
            label = "low_confidence"
        elif row["key_interaction_recall_consensus"] < config["triage"]["interaction_f1_threshold"]:
            label = "score_good_interactions_bad"
        else:
            label = "pending_review"
        labels.append(label)
    df["pose_decision_label"] = labels
    df["pose_decision_reason"] = df["pose_decision_label"]
    if "warnings_json" not in df.columns:
        df["warnings_json"] = df.get("warnings_json_x", "[]")
    out = df[["screening_pose_id", "molecule_id", "target_receptor_id", "receptor_state", "pose_decision_label", "pose_decision_reason", "pose_confidence_probability", "key_interaction_recall_consensus", "warnings_json"]]
    out.to_parquet(paths["processed"] / "screening_pose_decisions.parquet", index=False)
    out.to_csv(paths["processed"] / "screening_pose_decisions.csv", index=False)
    return out


def aggregate_candidates(conf: pd.DataFrame, decisions: pd.DataFrame, manifest: pd.DataFrame, master: pd.DataFrame, paths: dict) -> pd.DataFrame:
    df = conf.merge(decisions[["screening_pose_id", "pose_decision_label"]], on="screening_pose_id")
    man = manifest.drop_duplicates("molecule_id").set_index("molecule_id")
    master_idx = master.drop_duplicates("molecule_id").set_index("molecule_id")
    rows = []
    for mol, group in df.groupby("molecule_id"):
        best = group.sort_values("pose_confidence_probability", ascending=False).iloc[0]
        m = man.loc[mol]
        lib = master_idx.loc[mol]
        accepted_states = group[group["pose_decision_label"].eq("accepted_pose_consistent")]["receptor_state"].nunique()
        final_score = 0.55 * best["pose_confidence_probability"] + 0.25 * best["key_interaction_recall_consensus"] + 0.20 * (best["gnina_cnnscore"] if pd.notna(best["gnina_cnnscore"]) else 0)
        label = "accepted_candidate_control" if accepted_states else "low_confidence_or_interaction_rejected"
        rows.append({"molecule_id": mol, "source": m["source"], "subsource": m["subsource"], "screening_role": m["screening_role"], "screening_subset": m["screening_subset"], "standard_smiles": m["standard_smiles"], "novelty_bucket": m["novelty_bucket"], "closest_known_molecule_id": m["closest_known_molecule_id"], "tanimoto_to_closest_known": m["tanimoto_to_closest_known"], "scaffold_id": lib["scaffold_id"], "medchem_flags_json": m["medchem_flags_json"], "best_screening_pose_id": best["screening_pose_id"], "best_target_receptor_id": best["target_receptor_id"], "best_receptor_state": best["receptor_state"], "best_pose_confidence": best["pose_confidence_probability"], "best_calibrated_confidence": best["calibrated_confidence"], "best_docking_score": best["docking_score"], "best_gnina_cnnscore": best["gnina_cnnscore"], "best_gnina_cnnaffinity": best["gnina_cnnaffinity"], "best_ifp_tanimoto_to_consensus": best["ifp_tanimoto_to_consensus"], "best_key_interaction_recall_consensus": best["key_interaction_recall_consensus"], "num_receptors_screened": group["target_receptor_id"].nunique(), "num_receptor_states_with_accepted_pose": accepted_states, "receptor_state_consistency_score": accepted_states / max(group["receptor_state"].nunique(), 1), "source_normalized_rank": None, "novelty_adjusted_rank": None, "medchem_risk_score": 0.0 if m["medchem_flags_json"] == "[]" else 0.2, "final_candidate_score": final_score, "candidate_decision_label": label, "candidate_decision_reason": label, "warnings_json": json.dumps([])})
    out = pd.DataFrame(rows).sort_values("final_candidate_score", ascending=False)
    out["source_normalized_rank"] = out.groupby("source")["final_candidate_score"].rank(ascending=False, method="first")
    out["novelty_adjusted_rank"] = out.groupby("novelty_bucket")["final_candidate_score"].rank(ascending=False, method="first")
    out.to_parquet(paths["processed"] / "candidate_aggregate_scores.parquet", index=False)
    out.to_csv(paths["processed"] / "candidate_aggregate_scores.csv", index=False)
    return out
