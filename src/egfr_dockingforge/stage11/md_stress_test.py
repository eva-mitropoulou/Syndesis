from __future__ import annotations

from pathlib import Path
from typing import Any

import pandas as pd

from egfr_dockingforge.common.io import write_json
from egfr_dockingforge.stage11.candidate_selection import select_md_candidates
from egfr_dockingforge.stage11.forcefield_config import write_forcefield_policy
from egfr_dockingforge.stage11.gromacs_runner import make_gromacs_runs
from egfr_dockingforge.stage11.interaction_persistence import compute_interaction_persistence
from egfr_dockingforge.stage11.ligand_parameterization import parameterize_ligands
from egfr_dockingforge.stage11.load_stage_inputs import load_stage11_config, load_stage11_inputs, stage11_paths
from egfr_dockingforge.stage11.mdp_templates import write_mdp_templates
from egfr_dockingforge.stage11.report_stage11 import write_stage11_report
from egfr_dockingforge.stage11.stability_scoring import score_md_stability
from egfr_dockingforge.stage11.system_builder import build_md_systems
from egfr_dockingforge.stage11.trajectory_analysis import analyze_trajectories


def _load(config_path: str | Path) -> tuple[dict[str, Any], dict[str, Path], dict[str, pd.DataFrame | None]]:
    config = load_stage11_config(config_path)
    paths = stage11_paths(config)
    inputs = load_stage11_inputs(config)
    return config, paths, inputs


def select_md_candidates_cli(config_path: str | Path) -> dict:
    config, paths, inputs = _load(config_path)
    out = select_md_candidates(inputs, config, paths)
    write_forcefield_policy(config, paths)
    return {"status": "complete", "candidates": len(out)}


def parameterize_md_ligands_cli(config_path: str | Path) -> dict:
    config, paths, _ = _load(config_path)
    candidates = pd.read_parquet(paths["processed"] / "md_candidate_manifest.parquet")
    out = parameterize_ligands(candidates, config, paths)
    return {"status": "complete", "parameterized": int(out["parameterization_status"].eq("ready").sum()), "failed": int(out["parameterization_status"].ne("ready").sum())}


def build_md_systems_cli(config_path: str | Path) -> dict:
    config, paths, _ = _load(config_path)
    candidates = pd.read_parquet(paths["processed"] / "md_candidate_manifest.parquet")
    params = pd.read_parquet(paths["processed"] / "ligand_parameterization.parquet")
    write_mdp_templates(config, paths["md_root"] / "mdp_templates")
    out = build_md_systems(candidates, params, config, paths)
    return {"status": "complete", "systems": len(out), "ready": int(out["build_status"].isin(["ready_for_gromacs_build", "ready_for_md"]).sum())}


def run_md_minimization_cli(config_path: str | Path) -> dict:
    config, paths, _ = _load(config_path)
    systems = pd.read_parquet(paths["processed"] / "system_builds.parquet")
    runs = make_gromacs_runs(systems, config, paths)
    return {"status": "complete", "run_rows": len(runs)}


def run_md_equilibration_cli(config_path: str | Path) -> dict:
    return run_md_minimization_cli(config_path)


def run_md_production_cli(config_path: str | Path) -> dict:
    return run_md_minimization_cli(config_path)


def analyze_md_trajectories_cli(config_path: str | Path) -> dict:
    config, paths, _ = _load(config_path)
    runs = pd.read_parquet(paths["processed"] / "md_runs.parquet")
    metrics, qc = analyze_trajectories(runs, paths, config)
    return {"status": "complete", "metrics": len(metrics), "qc": len(qc)}


def compute_md_interaction_persistence_cli(config_path: str | Path) -> dict:
    config, paths, inputs = _load(config_path)
    candidates = pd.read_parquet(paths["processed"] / "md_candidate_manifest.parquet")
    metrics = pd.read_parquet(paths["processed"] / "md_metrics.parquet")
    long, summary = compute_interaction_persistence(candidates, metrics, inputs.get("key_interactions"), paths, config)
    return {"status": "complete", "interactions": len(long), "summary": len(summary)}


def score_md_stability_cli(config_path: str | Path) -> dict:
    config, paths, _ = _load(config_path)
    candidates = pd.read_parquet(paths["processed"] / "md_candidate_manifest.parquet")
    metrics = pd.read_parquet(paths["processed"] / "md_metrics.parquet")
    persistence = pd.read_parquet(paths["processed"] / "md_binding_mode_persistence.parquet")
    labels, summary, post = score_md_stability(candidates, metrics, persistence, paths)
    return {"status": "complete", "labels": len(labels), "post_md_rows": len(post)}


def report_stage11_cli(config_path: str | Path) -> dict:
    config, paths, _ = _load(config_path)
    report = write_stage11_report(paths)
    return {"status": "complete", "report": str(report)}


def run_stage11_all(config_path: str | Path) -> dict:
    select_md_candidates_cli(config_path)
    parameterize_md_ligands_cli(config_path)
    build_md_systems_cli(config_path)
    run_md_minimization_cli(config_path)
    analyze_md_trajectories_cli(config_path)
    compute_md_interaction_persistence_cli(config_path)
    score_md_stability_cli(config_path)
    summary = report_stage11_cli(config_path)
    config, paths, _ = _load(config_path)
    write_json(paths["processed"] / "stage11_summary.json", summary)
    return summary
