from __future__ import annotations

import html
from pathlib import Path

import pandas as pd


def _table(frame: pd.DataFrame, columns: list[str], n: int = 20) -> str:
    if frame.empty:
        return "<p>No rows.</p>"
    return frame[columns].head(n).to_html(index=False, escape=True)


def write_stage9_report(paths: dict[str, Path]) -> Path:
    seeds = pd.read_parquet(paths["processed"] / "analog_seed_scaffolds.parquet")
    edits = pd.read_parquet(paths["processed"] / "edit_sites.parquet")
    candidates = pd.read_parquet(paths["processed"] / "analog_candidates.parquet")
    validation = pd.read_parquet(paths["processed"] / "analog_validation.parquet")
    screening = pd.read_parquet(paths["processed"] / "analog_screening_results.parquet")
    acceptance = pd.read_parquet(paths["processed"] / "analog_acceptance.parquet")
    benchmark = pd.read_parquet(paths["processed"] / "analog_strategy_benchmark.parquet")
    provider = pd.read_parquet(paths["processed"] / "agent_proposal_status.parquet")
    target = paths["reports"] / "09_agentic_analog_optimization.html"
    target.write_text(
        f"""<!doctype html>
<html><head><meta charset=\"utf-8\"><title>Stage 9 Agentic Analog Optimization</title>
<style>body{{font-family:Arial,sans-serif;margin:32px;line-height:1.4}}table{{border-collapse:collapse;margin:16px 0}}td,th{{border:1px solid #bbb;padding:4px 7px}}th{{background:#eee}}</style></head>
<body>
<h1>Stage 9 - Interaction-constrained agentic analog optimization</h1>
<p><strong>Non-claim:</strong> all rows are computational proposals only. No experimental EGFR activity is claimed.</p>
<p>Agents propose transformations; RDKit, Stage 8 docking/GNINA/ProLIF, and Stage 6 pose confidence decide survival. Free-form LLM molecules are not accepted.</p>
<h2>Seed scaffolds</h2>{_table(seeds, ["seed_id","molecule_id","source","parent_candidate_rank","best_pose_confidence","best_key_interaction_recall_consensus"])}
<h2>Editable-site analysis</h2><p>{len(edits)} editable/protected site rows.</p>{_table(edits, ["edit_site_id","seed_id","attachment_atom_idx","protected_region_flag","editable_region_type"])}
<h2>Local LLM configuration and logs</h2><p>{html.escape(str(provider["proposal_status"].value_counts().to_dict()))}</p>{_table(provider, ["strategy_name","agent_role","proposal_status","warnings_json"])}
<h2>Analog generation and validation</h2><p>{len(candidates)} deterministic analog rows; {int(validation["hard_scope_pass"].sum())} passed hard scope.</p>{_table(validation, ["analog_id","hard_scope_pass","mw","clogp","tpsa","medchem_risk_score","rejection_reason"])}
<h2>Stage 8 mini-screen analog results</h2>{_table(screening, ["analog_id","strategy_name","best_pose_confidence","best_gnina_cnnscore","best_key_interaction_recall_consensus","binding_mode_preserved_flag"])}
<h2>Acceptance/rejection</h2>{_table(acceptance, ["analog_id","accepted_flag","acceptance_tier","delta_candidate_score","delta_pose_confidence","delta_ligand_efficiency","rejection_reason"])}
<h2>Benchmark</h2>{_table(benchmark, ["strategy_name","num_raw_proposals","num_valid_molecules","num_screened","num_accepted","accepted_analog_rate","benchmark_notes"])}
<h2>Failure taxonomy</h2><p>{html.escape(str(acceptance["acceptance_tier"].value_counts().to_dict()))}</p>
</body></html>
""",
        encoding="utf-8",
    )
    return target
