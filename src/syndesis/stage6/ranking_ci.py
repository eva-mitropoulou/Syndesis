"""Bootstrap confidence intervals for the pose-confidence ranking metrics.

The stage-6 headline (NDCG@3 0.41 -> 0.69) is a point estimate over a scaffold-holdout
TEST set with only ~10 ranking groups, so a reviewer rightly wants uncertainty. This
module recomputes NDCG@1/@3 and MRR on the TEST split for each model and attaches a
percentile bootstrap 95% CI by resampling the ranking GROUPS (the independent unit),
not individual poses. It reuses the exact per-group sklearn ndcg_score convention from
model_selection so the point estimates reproduce the reported numbers.

It also emits a "single best feature" honest-baseline comparison: the best individual
input feature (e.g. gnina_cnnscore) ranked directly, so the learned ranker's lift is
shown against the strongest trivial baseline, not only against raw docking score.
"""
from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from sklearn.metrics import ndcg_score


def _group_ndcg(group_scores: pd.DataFrame, k: int) -> float | None:
    """NDCG@k for one ranking group; None if the group has a single label value
    (undefined) -- matches model_selection which skips such groups."""
    if group_scores["label"].nunique() <= 1:
        return None
    y_true = group_scores["label"].to_numpy(dtype=float).reshape(1, -1)
    y_score = group_scores["score"].to_numpy(dtype=float).reshape(1, -1)
    return float(ndcg_score(y_true, y_score, k=min(k, len(group_scores))))


def _mrr(group_scores: pd.DataFrame) -> float | None:
    """Reciprocal rank of the first relevant (label>0) pose when ranked by score."""
    if (group_scores["label"] > 0).sum() == 0:
        return None
    ranked = group_scores.sort_values("score", ascending=False).reset_index(drop=True)
    for i, lab in enumerate(ranked["label"].to_numpy(), start=1):
        if lab > 0:
            return 1.0 / i
    return None


def _metric_over_groups(per_group: dict[str, pd.DataFrame], group_ids: list[str], metric: str) -> float:
    vals = []
    for gid in group_ids:
        g = per_group[gid]
        if metric == "NDCG@1":
            v = _group_ndcg(g, 1)
        elif metric == "NDCG@3":
            v = _group_ndcg(g, 3)
        elif metric == "MRR":
            v = _mrr(g)
        else:
            v = None
        if v is not None:
            vals.append(v)
    return float(np.mean(vals)) if vals else float("nan")


def bootstrap_ranking_ci(
    scored: pd.DataFrame,
    *,
    metrics: tuple[str, ...] = ("NDCG@1", "NDCG@3", "MRR"),
    n_boot: int = 3000,
    seed: int = 807,
) -> dict[str, Any]:
    """scored: rows with columns group_id, label, score (score = ranking signal,
    higher = better). Returns {metric: {point, ci_lo, ci_hi, n_groups}} with a
    percentile bootstrap over ranking groups."""
    per_group = {gid: g.reset_index(drop=True) for gid, g in scored.groupby("group_id")}
    group_ids = list(per_group.keys())
    rng = np.random.default_rng(seed)
    out: dict[str, Any] = {"n_groups": len(group_ids)}
    for metric in metrics:
        point = _metric_over_groups(per_group, group_ids, metric)
        boots = []
        n = len(group_ids)
        for _ in range(n_boot):
            sample = [group_ids[i] for i in rng.integers(0, n, n)]
            v = _metric_over_groups(per_group, sample, metric)
            if not np.isnan(v):
                boots.append(v)
        lo = float(np.percentile(boots, 2.5)) if boots else float("nan")
        hi = float(np.percentile(boots, 97.5)) if boots else float("nan")
        out[metric] = {
            "point": round(point, 4) if not np.isnan(point) else None,
            "ci_lo": round(lo, 4) if not np.isnan(lo) else None,
            "ci_hi": round(hi, 4) if not np.isnan(hi) else None,
        }
    return out


def build_test_scored_tables(paths: dict[str, Path], models: dict[str, str | Path]) -> dict[str, pd.DataFrame]:
    """For each (model_id -> artifact path), predict on the TEST split and return a
    {model_id: DataFrame[group_id,label,score]} table using the shared score_poses.
    Also builds baseline arms: raw docking score and each single input feature."""
    from syndesis.stage6.predict import score_poses

    features = pd.read_parquet(paths["processed"] / "pose_model_features.parquet")
    labels = pd.read_parquet(paths["processed"] / "pose_model_labels.parquet")
    splits = pd.read_parquet(paths["processed"] / "model_splits.parquet")

    # TEST poses (scaffold_holdout test fold)
    test_ids = set(splits.loc[splits["train_valid_test"] == "test", "pose_id"])
    feat = features[features["pose_id"].isin(test_ids)].copy()
    lab = labels[["pose_id", "rank_relevance_label"]].rename(columns={"rank_relevance_label": "label"})
    feat = feat.merge(lab, on="pose_id", how="left")
    feat["group_id"] = feat["docking_task_id"]

    tables: dict[str, pd.DataFrame] = {}

    # learned + configured model artifacts
    for model_id, artifact in models.items():
        try:
            preds = score_poses(feat, artifact)
            merged = feat[["pose_id", "group_id", "label"]].merge(
                preds[["pose_id", "pose_rank_model_score"]], on="pose_id", how="left"
            )
            merged = merged.rename(columns={"pose_rank_model_score": "score"})
            tables[model_id] = merged.dropna(subset=["score"])
        except Exception as exc:  # noqa: BLE001
            tables[f"{model_id}__FAILED"] = pd.DataFrame({"error": [str(exc)]})

    # honest baselines: raw docking score + single best features (higher=better)
    baseline_feats = {
        "baseline_docking_score": ("original_docking_score", False),  # more negative better -> invert
        "baseline_gnina_cnnscore": ("gnina_cnnscore", True),
        "baseline_gnina_cnnaffinity": ("gnina_cnnaffinity", True),
        "baseline_ifp_tanimoto_to_consensus": ("ifp_tanimoto_to_consensus", True),
        "baseline_key_interaction_recall_native": ("key_interaction_recall_native", True),
    }
    for name, (col, higher_better) in baseline_feats.items():
        if col not in feat.columns:
            continue
        t = feat[["pose_id", "group_id", "label"]].copy()
        s = pd.to_numeric(feat[col], errors="coerce")
        t["score"] = s if higher_better else -s
        tables[name] = t.dropna(subset=["score"])

    return tables


def run_ranking_ci(paths: dict[str, Path], models: dict[str, str | Path], out_dir: Path,
                   *, n_boot: int = 3000, seed: int = 807) -> pd.DataFrame:
    tables = build_test_scored_tables(paths, models)
    rows = []
    for arm, tbl in tables.items():
        if "error" in tbl.columns:
            rows.append({"arm": arm, "error": tbl["error"].iloc[0]})
            continue
        ci = bootstrap_ranking_ci(tbl, n_boot=n_boot, seed=seed)
        row: dict[str, Any] = {"arm": arm, "n_groups": ci["n_groups"]}
        for m in ("NDCG@1", "NDCG@3", "MRR"):
            row[f"{m}"] = ci[m]["point"]
            row[f"{m}_ci_lo"] = ci[m]["ci_lo"]
            row[f"{m}_ci_hi"] = ci[m]["ci_hi"]
        rows.append(row)
    report = pd.DataFrame(rows)
    out_dir.mkdir(parents=True, exist_ok=True)
    report.to_parquet(out_dir / "ranking_metrics_ci.parquet", index=False)
    report.to_csv(out_dir / "ranking_metrics_ci.csv", index=False)
    return report
