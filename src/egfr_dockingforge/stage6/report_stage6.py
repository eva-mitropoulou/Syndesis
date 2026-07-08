from __future__ import annotations

import html
from pathlib import Path

import pandas as pd


def _table(df: pd.DataFrame, max_rows: int = 20) -> str:
    return df.head(max_rows).to_html(index=False, escape=True)


def write_stage6_report(
    features: pd.DataFrame,
    labels: pd.DataFrame,
    groups: pd.DataFrame,
    audit: pd.DataFrame,
    splits: pd.DataFrame,
    metrics: pd.DataFrame,
    selection: pd.DataFrame,
    feature_importance: pd.DataFrame,
    ablation: pd.DataFrame,
    report_path: Path,
) -> Path:
    label_dist = labels["rank_relevance_label"].value_counts().sort_index().rename_axis("rank_relevance_label").reset_index(name="count")
    split_dist = splits.groupby(["split_name", "train_valid_test"]).size().reset_index(name="poses")
    leak = audit.groupby(["action", "leakage_risk"]).size().reset_index(name="features")
    tournament = metrics[metrics["metric_name"].isin(["valid_NDCG@1", "valid_NDCG@3", "valid_PR-AUC", "valid_Brier"])].pivot_table(
        index=["model_id", "model_family"], columns="metric_name", values="metric_value", aggfunc="max"
    ).reset_index().sort_values("valid_NDCG@1", ascending=False, na_position="last")
    html_text = f"""<!doctype html>
<html><head><meta charset="utf-8"><title>Stage 6 Pose Reranking Confidence</title>
<style>body{{font-family:Arial,sans-serif;margin:32px;line-height:1.4}}table{{border-collapse:collapse;margin:16px 0}}td,th{{border:1px solid #ccc;padding:4px 8px}}th{{background:#eee}}</style></head>
<body>
<h1>Stage 6 Pose Reranking and Confidence Calibration</h1>
<p>Rows: {len(features)} poses; labels: {len(labels)}; groups: {groups['group_id'].nunique()}.</p>
<h2>Label Distribution</h2>{_table(label_dist)}
<h2>Ranking Groups</h2>{_table(groups)}
<h2>Leakage Audit</h2>{_table(leak)}
<h2>Splits</h2>{_table(split_dist)}
<h2>Model Tournament</h2>{_table(tournament, 50)}
<h2>Model Selection</h2>{_table(selection)}
<h2>Feature Importance</h2>{_table(feature_importance, 30)}
<h2>Feature Group Ablation</h2>{_table(ablation, 30)}
<p>{html.escape(str(report_path))}</p>
</body></html>"""
    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text(html_text, encoding="utf-8")
    return report_path
