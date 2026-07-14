from __future__ import annotations

import json
import pandas as pd


def percentile(series: pd.Series, higher_is_better: bool) -> pd.Series:
    return series.rank(pct=True, ascending=higher_is_better, method="average")


def primary_triage(poses: pd.DataFrame, config: dict, paths: dict) -> pd.DataFrame:
    df = poses.copy()
    df["docking_score_percentile_within_receptor"] = df.groupby("target_receptor_id")["docking_score"].transform(lambda s: percentile(pd.to_numeric(s), False))
    df["docking_score_percentile_within_source"] = df.groupby("source")["docking_score"].transform(lambda s: percentile(pd.to_numeric(s), False))
    df["pass_primary_docking_triage"] = df.groupby(["molecule_id", "target_receptor_id"])["pose_rank"].transform(lambda s: s <= int(config["triage"]["top_n_per_molecule_receptor"]))
    df["triage_reason"] = "top_n_per_molecule_receptor"
    out = df[["screening_pose_id", "molecule_id", "target_receptor_id", "receptor_state", "pose_rank", "docking_score", "docking_score_percentile_within_receptor", "docking_score_percentile_within_source", "pass_primary_docking_triage", "triage_reason"]]
    out.to_parquet(paths["processed"] / "primary_docking_triage.parquet", index=False)
    out.to_csv(paths["processed"] / "primary_docking_triage.csv", index=False)
    return out


def normalize_scores(conf: pd.DataFrame, paths: dict) -> pd.DataFrame:
    df = conf.copy()
    df["docking_score_percentile"] = df.groupby("target_receptor_id")["docking_score"].transform(lambda s: percentile(pd.to_numeric(s), False))
    df["cnnscore_percentile"] = df.groupby("target_receptor_id")["gnina_cnnscore"].transform(lambda s: percentile(pd.to_numeric(s), True))
    df["cnnaffinity_percentile"] = df.groupby("target_receptor_id")["gnina_cnnaffinity"].transform(lambda s: percentile(pd.to_numeric(s), True))
    df["pose_confidence_percentile"] = df.groupby("target_receptor_id")["pose_confidence_probability"].transform(lambda s: percentile(pd.to_numeric(s), True))
    df["interaction_score_percentile"] = df.groupby("target_receptor_id")["key_interaction_recall_consensus"].transform(lambda s: percentile(pd.to_numeric(s), True))
    df["source_normalized_rank"] = df.groupby("docking_engine")["pose_confidence_probability"].rank(ascending=False, method="first")
    df["receptor_normalized_rank"] = df.groupby("target_receptor_id")["pose_confidence_probability"].rank(ascending=False, method="first")
    df["warnings_json"] = json.dumps([])
    out = df[["screening_pose_id", "molecule_id", "target_receptor_id", "receptor_state", "docking_score_percentile", "cnnscore_percentile", "cnnaffinity_percentile", "pose_confidence_percentile", "interaction_score_percentile", "source_normalized_rank", "receptor_normalized_rank", "warnings_json"]]
    out.to_parquet(paths["processed"] / "normalized_screening_scores.parquet", index=False)
    out.to_csv(paths["processed"] / "normalized_screening_scores.csv", index=False)
    return out
