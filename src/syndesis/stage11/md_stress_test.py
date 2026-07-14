from __future__ import annotations

from pathlib import Path
from typing import Any

import pandas as pd

from syndesis.common.io import write_json
from syndesis.stage11.candidate_selection import select_md_candidates
from syndesis.stage11.forcefield_config import write_forcefield_policy
from syndesis.stage11.gromacs_runner import make_gromacs_runs
from syndesis.stage11.interaction_persistence import compute_interaction_persistence
from syndesis.stage11.ligand_parameterization import parameterize_ligands
from syndesis.stage11.load_stage_inputs import load_stage11_config, load_stage11_inputs, stage11_paths
from syndesis.stage11.mdp_templates import write_mdp_templates
from syndesis.stage11.report_stage11 import write_stage11_report
from syndesis.stage11.stability_scoring import score_md_stability
from syndesis.stage11.system_builder import build_md_systems
from syndesis.stage11.trajectory_analysis import analyze_trajectories


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


def _assert_md_runs_succeeded(systems: pd.DataFrame, runs: pd.DataFrame, config: dict[str, Any] | None = None) -> None:
    """Fail loudly if any system that was ready_for_md did not produce its FULL
    set of completed production replicates. Two silent-failure modes this guards:
      * zero replicates completed (e.g. every system dies in equilibration on a
        bad .mdp) -> make_gromacs_runs still returns run rows and the driver would
        log a FALSE "OK run_md" over an empty/partial cohort;
      * fewer replicates than planned (e.g. selected_for_replicate_md not
        propagated, so 1 ran instead of 3) -> looks complete but is under-powered.
    Raising makes the driver mark the step FAILED (and retry after a fix)."""
    if systems is None or systems.empty or "build_status" not in systems.columns:
        return
    ready = systems[systems["build_status"] == "ready_for_md"]
    if ready.empty:
        return
    md_cfg = (config or {}).get("md", {}) if config else {}
    fin_reps = max(1, int(md_cfg.get("finalist_replicates", 3)))
    quick_reps = max(1, int(md_cfg.get("quick_replicates", 1)))

    def _planned(row: dict) -> int:
        if "planned_replicates" in row and pd.notna(row.get("planned_replicates")) and int(row["planned_replicates"]) > 0:
            return int(row["planned_replicates"])
        if bool(row.get("selected_for_replicate_md", False)):
            return fin_reps
        return quick_reps

    prod = runs[(runs["md_phase"] == "production_quick") & (runs["run_status"] == "complete")] \
        if runs is not None and not runs.empty and {"md_phase", "run_status", "md_system_id"} <= set(runs.columns) \
        else runs.iloc[0:0] if runs is not None else None
    done_counts = prod.groupby("md_system_id")["replicate_id"].nunique().to_dict() if prod is not None and not prod.empty and "replicate_id" in prod.columns else {}

    problems = []
    for row in ready.to_dict("records"):
        sid = row["md_system_id"]
        want = _planned(row)
        have = int(done_counts.get(sid, 0))
        if have < want:
            why = f"{have}/{want} production replicates complete"
            if have == 0 and runs is not None and not runs.empty and "md_system_id" in runs.columns:
                bad = runs[(runs["md_system_id"] == sid) & (runs["run_status"] != "complete")]
                if not bad.empty:
                    r = bad.iloc[0]
                    why += f" (furthest failure: {r.get('md_phase','?')} {r.get('run_status','?')})"
            problems.append(f"  - {sid}: {why}")
    if problems:
        raise RuntimeError(
            f"MD production incomplete for {len(problems)}/{len(ready)} ready systems:\n"
            + "\n".join(problems)
            + "\nFix the cause (bad .mdp / equilibration / replicate planning) and re-run; completed phases are skipped on resume."
        )


def run_md_minimization_cli(config_path: str | Path) -> dict:
    config, paths, _ = _load(config_path)
    systems = pd.read_parquet(paths["processed"] / "system_builds.parquet")
    # The replicate-planning flags live in the candidate manifest, not in
    # system_builds (the builder does not carry them). Without them,
    # _num_replicates() cannot see selected_for_replicate_md and silently falls
    # back to quick_replicates=1 -- so finalists that should get 3 independent
    # replicates only ran 1. Merge the planning columns in before running.
    manifest = pd.read_parquet(paths["processed"] / "md_candidate_manifest.parquet")
    plan_cols = [c for c in ["selected_for_quick_md", "selected_for_replicate_md", "planned_replicates"] if c in manifest.columns]
    if plan_cols:
        systems = systems.merge(
            manifest[["md_candidate_id", *plan_cols]].drop_duplicates("md_candidate_id"),
            on="md_candidate_id",
            how="left",
        )
    runs = make_gromacs_runs(systems, config, paths)
    _assert_md_runs_succeeded(systems, runs, config)
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
    labels, summary, post = score_md_stability(candidates, metrics, persistence, paths, config)
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
