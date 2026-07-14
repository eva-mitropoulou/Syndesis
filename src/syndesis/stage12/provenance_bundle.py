from __future__ import annotations

import hashlib
import json
import platform
import shutil
import subprocess
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable

from syndesis.common.io import resolve_path, write_json
from syndesis.stage12.evidence_aggregation import load_selection_and_inputs


def _sha256(path: Path) -> str | None:
    if not path.exists() or not path.is_file():
        return None
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _git_commit() -> str | None:
    try:
        return subprocess.check_output(["git", "rev-parse", "HEAD"], text=True).strip()
    except Exception:
        return None


def _copy_existing(paths: Iterable[Path], target_dir: Path) -> list[str]:
    copied: list[str] = []
    target_dir.mkdir(parents=True, exist_ok=True)
    for path in paths:
        if path.exists() and path.is_file():
            target = target_dir / path.name
            shutil.copy2(path, target)
            copied.append(str(target))
    return copied


def build_provenance_bundle(config_path: str | Path) -> dict[str, str]:
    config, paths, _, _ = load_selection_and_inputs(config_path)
    bundle = paths["bundle"]
    config_files = [resolve_path(config_path)]
    source_registry = resolve_path(config["inputs"]["stage12_sources"])
    output_tables = [paths["processed"] / "final_candidate_selection.parquet", paths["processed"] / "final_ranked_candidates.parquet"]
    report_files = [
        paths["reports"] / "12_final_candidate_summary.html",
        paths["reports"] / "final_candidate_table.csv",
        paths["model_cards"] / "pose_confidence_model_card.md",
        paths["dataset_cards"] / "project_dataset_card.md",
    ]
    _copy_existing(config_files, bundle / "configs")
    _copy_existing([source_registry], bundle / "source_registries")
    _copy_existing(output_tables, bundle / "tables")
    _copy_existing(report_files, bundle / "reports")
    _copy_existing(paths["cards"].glob("*.json"), bundle / "candidate_cards")
    input_hashes = {
        name: _sha256(resolve_path(path))
        for name, path in config.get("inputs", {}).items()
        if name != "stage12_sources" and resolve_path(path).exists()
    }
    output_hashes = {str(path): _sha256(path) for path in output_tables if path.exists()}
    manifest = {
        "project_name": "syndesis",
        "project_version": "stage12.v1",
        "run_id": f"stage12_{datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')}",
        "git_commit": _git_commit(),
        "environment": {"python": platform.python_version(), "platform": platform.platform()},
        "config_hashes": {str(path): _sha256(path) for path in config_files},
        "input_table_hashes": input_hashes,
        "output_table_hashes": output_hashes,
        "model_artifacts": ["reports/model_cards/pose_confidence_model_card.md"],
        "report_files": [str(path) for path in report_files if path.exists()],
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "tool_versions": {"stage12": "python_cli", "rdkit": "optional_for_svg_depictions"},
        "command_log": [
            "syndesis build-final-candidate-table --config configs/stage12_candidate_dossiers.yaml",
            "syndesis build-candidate-cards --config configs/stage12_candidate_dossiers.yaml",
            "syndesis render-candidate-dossiers --config configs/stage12_candidate_dossiers.yaml",
            "syndesis render-final-figures --config configs/stage12_candidate_dossiers.yaml",
            "syndesis build-model-card --config configs/stage12_candidate_dossiers.yaml",
            "syndesis build-dataset-card --config configs/stage12_candidate_dossiers.yaml",
            "syndesis report-stage12 --config configs/stage12_candidate_dossiers.yaml",
        ],
        "limitations": [
            "No experimental activity measurements are included.",
            "Current final table separates diagnostic known controls from rejected generated analog examples.",
            "Stage 11 production MD trajectory evidence is not available in this run.",
        ],
    }
    manifest_path = bundle / "manifest.json"
    write_json(manifest_path, manifest)
    return {"manifest": str(manifest_path)}
