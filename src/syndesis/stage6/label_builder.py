from __future__ import annotations

import numpy as np
import pandas as pd

from syndesis.common.io import write_table
from syndesis.stage6 import schemas


RELEVANCE = {
    "high_confidence_native_like": 3,
    "plausible_binding_mode": 2,
    "rmsd_poor_interactions_good": 1,
    "uncertain_plausible_pose": 1,
    "wrong_binding_mode": 0,
    "invalid_pose": 0,
    "sampling_failure": 0,
    "failed_analysis": 0,
}


def build_pose_labels(inputs: dict[str, pd.DataFrame], config: dict, paths: dict) -> pd.DataFrame:
    final = inputs["final_pose_labels"].copy()
    stage3 = inputs["stage3_pose_labels"][["pose_id", "stage3_pose_label", "strict_native_like_flag", "relaxed_native_like_flag", "invalid_pose_flag"]].copy()
    labels = final.merge(stage3, on="pose_id", how="left", suffixes=("", "_stage3"))

    # Label definition. The pose-confidence model predicts whether a docked pose
    # reproduces the crystallographic (native) binding geometry. The DEFAULT and
    # scientifically-clean target is therefore the RMSD-to-native ground truth
    # (strict/relaxed native-like flags from Stage 3), NOT the Stage 5
    # interaction-consistency assessment. Using the interaction assessment as the
    # label while also feeding interaction-recovery features to the model creates
    # target leakage (the model reads a transform of its own label); those
    # features are now dropped by the leakage audit, and the label is decoupled
    # here so the reported skill reflects learning from pose-intrinsic features.
    label_definition = str(config["labels"].get("label_definition", "rmsd_native_like"))
    strict = labels["strict_native_like_flag"].fillna(False).astype(bool)
    relaxed = labels["relaxed_native_like_flag"].fillna(False).astype(bool)
    invalid = labels["invalid_pose_flag"].fillna(False).astype(bool)

    if label_definition == "rmsd_native_like":
        # 3 = strict native-like (tight RMSD), 2 = relaxed native-like, 0 = other.
        rank = np.where(strict, 3, np.where(relaxed, 2, 0))
        labels["rank_relevance_label"] = rank
        labels.loc[invalid, "rank_relevance_label"] = 0
        labels["rank_relevance_label"] = labels["rank_relevance_label"].astype(int)
        labels["label_reason"] = np.where(
            labels["rank_relevance_label"].ge(2),
            "native-like pose by RMSD-to-crystal ground truth (Stage 3)",
            "non-native-like or invalid pose by RMSD-to-crystal ground truth",
        )
    else:
        # Legacy hybrid definition (interaction-consistency + RMSD flags).
        labels["rank_relevance_label"] = labels["final_pose_label"].map(RELEVANCE).fillna(0).astype(int)
        labels.loc[strict, "rank_relevance_label"] = 3
        labels.loc[relaxed & labels["rank_relevance_label"].lt(2), "rank_relevance_label"] = 2
        labels.loc[invalid, "rank_relevance_label"] = 0
        labels["label_reason"] = np.where(
            labels["rank_relevance_label"].ge(1),
            "interaction-consistent pose label from Stage 5 (legacy hybrid)",
            "wrong/invalid/failed pose label from Stage 5 or Stage 3 (legacy hybrid)",
        )

    min_positive = int(config["labels"].get("binary_positive_min_relevance", 1))
    labels["binary_confidence_label"] = labels["rank_relevance_label"].ge(min_positive).astype(int)
    confidence_map = {"high": 1.0, "medium": 0.7, "low": 0.4}
    labels["label_confidence"] = labels["final_pose_label_confidence"].map(confidence_map).fillna(0.5)
    labels["stage3_pose_label"] = labels["stage3_pose_label"].fillna(labels.get("stage3_pose_label_stage3"))
    out = labels[schemas.LABEL_COLUMNS].copy()
    write_table(paths["processed"] / "pose_model_labels.parquet", out)
    write_table(paths["processed"] / "pose_model_labels.csv", out)
    return out
