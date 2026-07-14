from __future__ import annotations

from pathlib import Path
from typing import Any

import pandas as pd

from syndesis.common.io import ensure_dir, write_json
from syndesis.stage3.docking_prep import prepare_ligands, prepare_receptors
from syndesis.stage3.load_stage2 import load_stage2_ensemble, load_stage3_config, stage3_paths
from syndesis.stage3.pose_labeling import build_labels, task_metrics
from syndesis.stage3.pose_sanity import sanity_for_poses
from syndesis.stage3.receptor_validation import receptor_validation_summary
from syndesis.stage3.report_stage3 import write_stage3_report
from syndesis.stage3.rmsd import compute_pose_rmsd
from syndesis.stage3.task_matrix import build_reference_transforms, build_task_matrix
from syndesis.stage3.vina_runner import run_vina_tasks
from syndesis.stage3.unidock_runner import run_unidock_tasks


def ensure_stage3_dirs(paths: dict[str, Path]) -> None:
    for path in paths.values():
        ensure_dir(path)


def prepare_docking_inputs(config_path: str | Path) -> dict[str, Any]:
    config = load_stage3_config(config_path)
    paths = stage3_paths(config)
    ensure_stage3_dirs(paths)
    ensemble = load_stage2_ensemble(config)
    receptors = prepare_receptors(ensemble, config, paths)
    ligands = prepare_ligands(ensemble, config, paths)
    return {"status": "complete", "receptors": int(len(receptors)), "ligands": int(len(ligands))}


def build_docking_task_matrix(config_path: str | Path) -> dict[str, Any]:
    config = load_stage3_config(config_path)
    paths = stage3_paths(config)
    ensure_stage3_dirs(paths)
    ensemble = load_stage2_ensemble(config)
    receptor_prep = pd.read_parquet(paths["processed"] / "receptor_docking_prep.parquet") if (paths["processed"] / "receptor_docking_prep.parquet").exists() else prepare_receptors(ensemble, config, paths)
    ligand_prep = pd.read_parquet(paths["processed"] / "ligand_docking_prep.parquet") if (paths["processed"] / "ligand_docking_prep.parquet").exists() else prepare_ligands(ensemble, config, paths)
    transforms = build_reference_transforms(ensemble, ligand_prep, paths)
    tasks = build_task_matrix(ensemble, receptor_prep, ligand_prep, transforms, config, paths)
    return {"status": "complete", "tasks": int(len(tasks))}


def run_redocking(config_path: str | Path) -> dict[str, Any]:
    return run_all_docking(config_path)


def run_crossdocking(config_path: str | Path) -> dict[str, Any]:
    return run_all_docking(config_path)


def run_all_docking(config_path: str | Path) -> dict[str, Any]:
    config = load_stage3_config(config_path)
    paths = stage3_paths(config)
    task_path = paths["processed"] / "docking_task_matrix.parquet"
    if not task_path.exists():
        build_docking_task_matrix(config_path)
    tasks = pd.read_parquet(task_path)
    if config["docking"]["primary_engine"] == "unidock":
        runs, poses = run_unidock_tasks(tasks, config, paths)
    else:
        runs, poses = run_vina_tasks(tasks, config, paths)
    return {"status": "complete", "runs": int(len(runs)), "poses": int(len(poses))}


def compute_docking_rmsd(config_path: str | Path) -> dict[str, Any]:
    config = load_stage3_config(config_path)
    paths = stage3_paths(config)
    poses = pd.read_parquet(paths["processed"] / "docked_poses.parquet")
    tasks = pd.read_parquet(paths["processed"] / "docking_task_matrix.parquet")
    transforms = pd.read_parquet(paths["processed"] / "reference_pose_transforms.parquet")
    ligand_prep_path = paths["processed"] / "ligand_docking_prep.parquet"
    ligand_prep = pd.read_parquet(ligand_prep_path) if ligand_prep_path.exists() else None
    rmsd = compute_pose_rmsd(poses, tasks, transforms, config, paths, ligand_prep)
    return {"status": "complete", "poses": int(len(rmsd))}


def run_pose_sanity_checks(config_path: str | Path) -> dict[str, Any]:
    config = load_stage3_config(config_path)
    paths = stage3_paths(config)
    poses = pd.read_parquet(paths["processed"] / "docked_poses.parquet")
    sanity = sanity_for_poses(poses, paths)
    return {"status": "complete", "poses": int(len(sanity))}


def label_stage3_poses(config_path: str | Path) -> dict[str, Any]:
    config = load_stage3_config(config_path)
    paths = stage3_paths(config)
    rmsd = pd.read_parquet(paths["processed"] / "pose_rmsd.parquet")
    sanity = pd.read_parquet(paths["processed"] / "pose_sanity.parquet")
    labels = build_labels(rmsd, sanity, config, paths)
    tasks = pd.read_parquet(paths["processed"] / "docking_task_matrix.parquet")
    runs = pd.read_parquet(paths["processed"] / "docking_runs.parquet")
    metrics = task_metrics(tasks, labels, runs, config, paths)
    ensemble = load_stage2_ensemble(config)
    validation = receptor_validation_summary(ensemble, metrics, paths)
    return {"status": "complete", "labels": int(len(labels)), "task_metrics": int(len(metrics)), "receptor_validation": int(len(validation))}


def report_stage3(config_path: str | Path) -> dict[str, Any]:
    config = load_stage3_config(config_path)
    paths = stage3_paths(config)
    tasks = pd.read_parquet(paths["processed"] / "docking_task_matrix.parquet")
    runs = pd.read_parquet(paths["processed"] / "docking_runs.parquet")
    metrics = pd.read_parquet(paths["processed"] / "docking_task_metrics.parquet")
    validation = pd.read_parquet(paths["processed"] / "receptor_stage3_validation.parquet")
    report = write_stage3_report(tasks, runs, metrics, validation, paths["reports"] / "03_redocking_crossdocking.html")
    summary = {"status": "complete", "tasks": int(len(tasks)), "runs": int(len(runs)), "report": str(report)}
    write_json(paths["processed"] / "stage3_summary.json", summary)
    return summary


def build_stage3_all(config_path: str | Path) -> dict[str, Any]:
    prepare_docking_inputs(config_path)
    build_docking_task_matrix(config_path)
    run_all_docking(config_path)
    compute_docking_rmsd(config_path)
    run_pose_sanity_checks(config_path)
    label_stage3_poses(config_path)
    return report_stage3(config_path)
