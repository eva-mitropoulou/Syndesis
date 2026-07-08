from __future__ import annotations

import json

import joblib
import pandas as pd

from egfr_dockingforge.stage6.predict import score_poses


def apply_pose_confidence(poses: pd.DataFrame, gnina: pd.DataFrame, sanity: pd.DataFrame, interactions: pd.DataFrame, tasks: pd.DataFrame, manifest: pd.DataFrame, config: dict, paths: dict) -> tuple[pd.DataFrame, pd.DataFrame]:
    task_meta = tasks.set_index("screening_task_id")
    rows = []
    merged = poses.merge(gnina, on=["screening_pose_id", "molecule_id", "target_receptor_id", "receptor_state", "docking_score"], how="inner")
    merged = merged.merge(sanity[["screening_pose_id", "sanity_status"]], on="screening_pose_id", how="left")
    merged = merged.merge(interactions, on=["screening_pose_id", "molecule_id", "target_receptor_id", "receptor_state"], how="left")
    man = manifest.set_index("prepared_ligand_id")
    for row in merged.to_dict("records"):
        task = task_meta.loc[row["screening_task_id"]]
        mrow = man.loc[row["prepared_ligand_id"]]
        rows.append({"pose_id": row["screening_pose_id"], "screening_pose_id": row["screening_pose_id"], "docking_task_id": row["screening_task_id"], "ligand_id": row["molecule_id"], "target_receptor_id": row["target_receptor_id"], "task_type": "stage8_screening", "receptor_state": row["receptor_state"], "docking_engine": row["docking_engine"], "original_pose_rank": row["pose_rank"], "pose_rank": row["pose_rank"], "original_docking_score": row["docking_score"], "docking_score": row["docking_score"], "cnnscore": row["cnnscore"], "cnnaffinity": row["cnnaffinity"], "gnina_cnnscore": row["cnnscore"], "gnina_cnnaffinity": row["cnnaffinity"], "gnina_empirical_affinity": row["gnina_empirical_affinity"], "cnn_vs": row["cnn_vs"], "sanity_status": row["sanity_status"], "ifp_tanimoto_to_consensus": row["ifp_tanimoto_to_consensus"], "key_interaction_recall_consensus": row["key_interaction_recall_consensus"], "key_interaction_precision_consensus": row["key_interaction_precision_consensus"], "key_interaction_f1_consensus": row["key_interaction_f1_consensus"], "num_interactions": row["num_interactions"], "num_key_interactions": row["num_key_interactions"], "fingerprint_sparse_json": row["fingerprint_sparse_json"], "binding_mode_cluster_id": row["binding_mode_cluster_id"], "ligand_source": mrow["source"], "receptor_cluster_id": task["receptor_cluster_id"], "rank_fraction_within_task": row["pose_rank"] / max(int(task["num_modes"]), 1), "sanity_status_encoded": 1.0 if row["sanity_status"] == "pass" else 0.0, "binding_mode_cluster_distance": 1.0 - row["binding_mode_compatibility_score"], "dominant_binding_mode_cluster_compatibility": row["binding_mode_compatibility_score"], "state_match_flag": 1, "docking_box_volume": float(task["docking_box_size_x"]) * float(task["docking_box_size_y"]) * float(task["docking_box_size_z"])})
    features = pd.DataFrame(rows)
    artifact = joblib.load(config["inputs"]["stage6_model"])
    scored = score_poses(features, artifact)
    conf = features.merge(scored, left_on="screening_pose_id", right_on="pose_id", how="left")
    out = pd.DataFrame({"screening_pose_id": conf["screening_pose_id"], "molecule_id": conf["ligand_id"], "target_receptor_id": conf["target_receptor_id"], "receptor_state": conf["receptor_state"], "docking_engine": conf["docking_engine"], "original_pose_rank": conf["original_pose_rank"], "docking_score": conf["docking_score"], "gnina_cnnscore": conf["gnina_cnnscore"], "gnina_cnnaffinity": conf["gnina_cnnaffinity"], "sanity_status": conf["sanity_status"], "ifp_tanimoto_to_consensus": conf["ifp_tanimoto_to_consensus"], "key_interaction_recall_consensus": conf["key_interaction_recall_consensus"], "binding_mode_cluster_id": conf["binding_mode_cluster_id"], "pose_rank_model_score": conf["pose_rank_model_score"], "pose_confidence_probability": conf["pose_confidence_probability"], "model_rank_within_task": conf["model_rank_within_group"], "calibrated_confidence": conf["pose_confidence_probability"], "confidence_label": conf["confidence_label"], "model_id": artifact["selected_ranker_model_id"], "model_version": "stage6_pose_ranker_confidence_v1", "prediction_status": "success", "warnings_json": json.dumps([])})
    features.to_parquet(paths["processed"] / "screening_pose_features.parquet", index=False)
    features.to_csv(paths["processed"] / "screening_pose_features.csv", index=False)
    out.to_parquet(paths["processed"] / "screening_pose_confidence.parquet", index=False)
    out.to_csv(paths["processed"] / "screening_pose_confidence.csv", index=False)
    return features, out
