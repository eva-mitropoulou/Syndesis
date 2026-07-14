from __future__ import annotations

from pathlib import Path

import pandas as pd


def _table(df: pd.DataFrame, cols: list[str], n: int = 30) -> str:
    if df.empty:
        return "<p>No rows.</p>"
    return df[cols].head(n).to_html(index=False, escape=True)


def write_stage11_report(paths: dict[str, Path]) -> Path:
    candidates = pd.read_parquet(paths["processed"] / "md_candidate_manifest.parquet")
    params = pd.read_parquet(paths["processed"] / "ligand_parameterization.parquet")
    param_report = pd.read_parquet(paths["processed"] / "ligand_parameterization_report.parquet")
    systems = pd.read_parquet(paths["processed"] / "system_builds.parquet")
    runs = pd.read_parquet(paths["processed"] / "md_runs.parquet")
    labels = pd.read_parquet(paths["processed"] / "md_pose_stability_labels.parquet")
    summary = pd.read_parquet(paths["processed"] / "md_candidate_summary.parquet")
    target = paths["reports"] / "11_md_stress_test.html"
    target.write_text(
        f"""<!doctype html>
<html><head><meta charset="utf-8"><title>Stage 11 MD Stress Test</title>
<style>body{{font-family:Arial,sans-serif;margin:32px;line-height:1.4}}table{{border-collapse:collapse;margin:16px 0}}td,th{{border:1px solid #bbb;padding:4px 7px}}th{{background:#eee}}</style></head>
<body>
<h1>Stage 11 - Explicit-solvent all-atom MD stress testing</h1>
<p><strong>Non-claim:</strong> MD is a computational stress test and is not experimental proof of EGFR binding, potency, selectivity, or clinical relevance.</p>
<p>The configured primary stack is GROMACS + AMBER protein force field + GAFF2 ligand parameters generated with AmberTools/ACPYPE and AM1-BCC charges. GAFF2 parameters are an open, practical workflow, not a guarantee of perfect ligand physics; parameterization warnings must be reviewed.</p>
<h2>Candidate selection</h2>{_table(candidates, ["md_candidate_id","molecule_id","source","stage8_candidate_score","selected_for_quick_md","selected_for_replicate_md","selection_reason"])}
<h2>Parameterization</h2>{_table(params, ["md_candidate_id","parameterization_status","rejection_reason","cgenff_penalty_status","charmm2gmx_status","warnings_json"])}
<h2>AMBER/GAFF2 parameterization report</h2>{_table(param_report, ["ligand_id","backend","net_charge","charge_model","ligand_forcefield","protein_forcefield","water_model","parameterization_status","warnings_json"])}
<h2>System builds</h2>{_table(systems, ["md_system_id","md_candidate_id","water_model","box_type","salt_concentration_molar","build_status","warnings_json"])}
<h2>MD runs</h2>{_table(runs, ["md_run_id","md_phase","ensemble","gromacs_version","run_status","error_message"])}
<h2>Stability labels</h2>{_table(labels, ["md_candidate_id","replicate_id","md_stability_label","md_stability_score","md_rejection_reason"])}
<h2>Candidate summary</h2>{_table(summary, ["md_candidate_id","final_md_label","final_md_decision","final_md_reason","recommended_for_final_dossier"])}
<h2>Limitations</h2><p>Ligand parameters use an open AMBER/GAFF2/AM1-BCC workflow and should be inspected before flagship MD. Short-timescale MD remains a stress test rather than a binding free-energy or experimental activity claim.</p>
</body></html>
""",
        encoding="utf-8",
    )
    return target
