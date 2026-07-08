from __future__ import annotations

import html
from pathlib import Path

import pandas as pd


def _table(frame: pd.DataFrame, columns: list[str], n: int = 30) -> str:
    if frame.empty:
        return "<p>No rows.</p>"
    return frame[columns].head(n).to_html(index=False, escape=True)


def write_stage10_report(paths: dict[str, Path]) -> Path:
    manifest = pd.read_parquet(paths["processed"] / "ablation_strategy_manifest.parquet")
    budget = pd.read_parquet(paths["processed"] / "strategy_budget_audit.parquet")
    metrics = pd.read_parquet(paths["processed"] / "strategy_metrics.parquet")
    stats = pd.read_parquet(paths["processed"] / "statistical_comparisons.parquet")
    ablations = pd.read_parquet(paths["processed"] / "ablation_summary.parquet")
    hacking = pd.read_parquet(paths["processed"] / "score_hacking_cases.parquet")
    diversity = pd.read_parquet(paths["processed"] / "diversity_novelty_metrics.parquet")
    cost = pd.read_parquet(paths["processed"] / "compute_cost_metrics.parquet")
    target = paths["reports"] / "10_ablation_benchmark.html"
    target.write_text(
        f"""<!doctype html>
<html><head><meta charset="utf-8"><title>Stage 10 Ablation Benchmark</title>
<style>body{{font-family:Arial,sans-serif;margin:32px;line-height:1.4}}table{{border-collapse:collapse;margin:16px 0}}td,th{{border:1px solid #bbb;padding:4px 7px}}th{{background:#eee}}</style></head>
<body>
<h1>Stage 10 - Budget-controlled analog optimization benchmark</h1>
<p><strong>Non-claim:</strong> this report is a computational benchmark only and makes no experimental EGFR activity claim.</p>
<p>Stage 10 consumed Stage 9 outputs only. No new docking, GNINA, ProLIF, pose-confidence scoring, agent generation, or MD was run.</p>
<h2>Benchmark design</h2><p>Paired seed scaffolds, fixed proposal/screening budgets, accepted analog rate as the primary metric, score-hacking rate as the primary failure metric.</p>
<h2>Strategy definitions</h2>{_table(manifest, ["strategy_id","strategy_name","enabled_flag","uses_llm","uses_prolif_constraint","uses_pose_confidence","uses_md_filter","md_status"])}
<h2>Budget audit</h2>{_table(budget, ["strategy_id","seed_id","num_raw_proposals","num_valid_unique_analogs","num_screened_analogs","budget_normalized_flag","budget_violation_reason"])}
<h2>Pre-MD benchmark results</h2>{_table(metrics, ["strategy_name","num_raw_proposals","num_valid_molecules","num_screened","num_pre_md_accepted","accepted_analog_rate_pre_md","score_hacking_rate","md_status"])}
<h2>Post-MD benchmark status</h2><p>{html.escape(str(manifest["md_status"].value_counts().to_dict()))}</p>
<h2>Ablation summary</h2>{_table(ablations, ["ablation_id","base_strategy","added_component","accepted_rate_change","score_hacking_rate_change","conclusion","evidence_strength"])}
<h2>Score-hacking analysis</h2>{_table(hacking, ["analog_id","strategy_id","score_hacking_type","improved_metric","worsened_metric","severity"])}
<h2>Statistical comparisons</h2>{_table(stats, ["comparison_id","delta_mean","ci_low","ci_high","p_value_corrected","interpretation"])}
<h2>Diversity and novelty</h2>{_table(diversity, ["strategy_id","seed_id","internal_diversity","mean_parent_tanimoto","mode_collapse_flag"])}
<h2>Compute and cost</h2>{_table(cost, ["strategy_id","num_docking_tasks","num_gnina_tasks","num_prolif_tasks","accepted_analogs_per_gpu_hour","notes"])}
<h2>Limitations</h2><p>Stage 10 is a computational benchmark over the available Stage 9 strategy outputs. Statistical outputs are exploratory and should not be interpreted as experimental significance claims.</p>
</body></html>
""",
        encoding="utf-8",
    )
    return target
