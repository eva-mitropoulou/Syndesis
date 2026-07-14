from __future__ import annotations

import json
import subprocess
import time
from pathlib import Path
from typing import Any

import pandas as pd

from syndesis.common.io import ensure_dir, project_root, write_table
from syndesis.stage4.score_parser import parse_gnina_output
from syndesis.stage4.scoring_engines import check_gnina
from syndesis.stage4.schemas import GNINA_SCORE_COLUMNS, RAW_GNINA_COLUMNS


def _relativize(path: str | Path) -> str:
    root = project_root()
    candidate = Path(path)
    try:
        return str(candidate.relative_to(root))
    except ValueError:
        return str(candidate)


def _docker_gnina_engine(gnina_cfg: dict[str, Any]):
    """Probe the gnina/gnina Docker image and return an EngineInfo-like record."""
    from syndesis.stage4.scoring_engines import EngineInfo, parse_gnina_version

    image = str(gnina_cfg.get("docker_image", "gnina/gnina:latest"))
    try:
        completed = subprocess.run(
            ["docker", "run", "--rm", "--gpus", "all", image, "gnina", "--version"],
            check=False, capture_output=True, text=True, timeout=120,
        )
    except Exception as exc:  # noqa: BLE001
        return EngineInfo(False, f"docker:{image}", None, f"gnina docker probe failed: {exc}")
    text = f"{completed.stdout}\n{completed.stderr}"
    version = parse_gnina_version(text)
    return EngineInfo(completed.returncode == 0, f"docker:{image}", version, None if completed.returncode == 0 else text[-1000:])


def build_gnina_command(gnina_cfg: dict[str, Any], gnina_args: list[str]) -> list[str]:
    """Build a GNINA command, optionally routed through the gnina/gnina Docker
    image. When ``use_docker`` is true the project root is bind-mounted at /work
    so the relative -r/-l/-o paths resolve inside the container. All GNINA file
    arguments must therefore be project-root-relative (use ``_relativize``)."""
    if gnina_cfg.get("use_docker"):
        image = str(gnina_cfg.get("docker_image", "gnina/gnina:latest"))
        root = str(project_root())
        return [
            "docker", "run", "--rm", "--gpus", "all",
            "-v", f"{root}:/work", "-w", "/work",
            image, "gnina", *gnina_args,
        ]
    return [str(gnina_cfg.get("executable", "gnina")), *gnina_args]


def _prepare_ligand_for_gnina(pose_file: str | Path, paths: dict[str, Path], obabel: str | None = None) -> Path:
    """Convert a docked pose to an SDF (with bond orders) for GNINA scoring.

    GNINA's CNN and empirical scoring need correct connectivity/bond orders.
    A PDBQT stripped to bare ATOM/HETATM lines loses all bond information, so
    the previous PDB line-strip corrupted scoring. We use OpenBabel to convert
    PDBQT -> SDF, which perceives and writes bonds. Conversion failure blocks
    rescoring because a bond-order-deficient substitute is not scientifically
    equivalent.
    """
    source = Path(pose_file)
    if source.suffix.lower() not in {".pdbqt", ".pdb"}:
        return source
    ligand_dir = ensure_dir(paths["processed"] / "gnina_ligands")
    target = ligand_dir / f"{source.stem}.sdf"
    if target.exists() and target.stat().st_mtime >= source.stat().st_mtime:
        return target
    if not obabel:
        raise RuntimeError("Open Babel is required to convert PDB/PDBQT poses for GNINA rescoring.")
    completed = subprocess.run(
        [obabel, str(source), "-O", str(target)],
        check=False,
        capture_output=True,
        text=True,
        timeout=120,
    )
    if completed.returncode != 0 or not target.exists() or target.stat().st_size == 0:
        detail = (completed.stderr or completed.stdout or "no conversion output").strip()
        raise RuntimeError(f"Open Babel ligand conversion failed for {source}: {detail[-500:]}")
    return target


def _rank_within_task(frame: pd.DataFrame, column: str, ascending: bool) -> pd.Series:
    values = pd.to_numeric(frame[column], errors="coerce")
    return values.groupby(frame["docking_task_id"]).rank(method="min", ascending=ascending, na_option="bottom")


def _write_gnina_outputs(raw: pd.DataFrame, scores: pd.DataFrame, paths: dict[str, Path]) -> None:
    if not scores.empty:
        scores["cnnscore_rank_within_task"] = _rank_within_task(scores, "cnnscore", ascending=False)
        scores["cnnaffinity_rank_within_task"] = _rank_within_task(scores, "cnnaffinity", ascending=False)
        scores["gnina_affinity_rank_within_task"] = _rank_within_task(scores, "gnina_empirical_affinity", ascending=True)
    write_table(paths["processed"] / "gnina_raw_runs.parquet", raw)
    write_table(paths["processed"] / "gnina_raw_runs.csv", raw)
    write_table(paths["processed"] / "gnina_scores.parquet", scores)
    write_table(paths["processed"] / "gnina_scores.csv", scores)


def run_gnina_rescoring(tasks: pd.DataFrame, config: dict[str, Any], paths: dict[str, Path]) -> tuple[pd.DataFrame, pd.DataFrame]:
    gnina_cfg = config["gnina"]
    if gnina_cfg.get("use_docker"):
        engine = _docker_gnina_engine(gnina_cfg)
    else:
        engine = check_gnina(str(gnina_cfg.get("executable", "gnina")))
    ensure_dir(paths["logs"])
    raw_path = paths["processed"] / "gnina_raw_runs.parquet"
    scores_path = paths["processed"] / "gnina_scores.parquet"
    raw = pd.read_parquet(raw_path) if raw_path.exists() else pd.DataFrame(columns=RAW_GNINA_COLUMNS)
    scores = pd.read_parquet(scores_path) if scores_path.exists() else pd.DataFrame(columns=GNINA_SCORE_COLUMNS)
    completed = set(scores.loc[scores["rescoring_status"].eq("success"), "rescoring_task_id"]) if not scores.empty else set()
    raw_rows: list[dict[str, Any]] = raw.to_dict("records")
    score_rows: list[dict[str, Any]] = scores.to_dict("records")

    ready = tasks[tasks["task_status"] == "ready"].copy()
    for _, task in ready.iterrows():
        task_id = str(task["rescoring_task_id"])
        if task_id in completed:
            continue
        stdout_log = paths["logs"] / f"{task_id}.stdout.log"
        stderr_log = paths["logs"] / f"{task_id}.stderr.log"
        ligand_file = _prepare_ligand_for_gnina(task["pose_file"], paths, gnina_cfg.get("obabel") or config.get("prep", {}).get("obabel_path"))
        gnina_args = ["--score_only", "-r", _relativize(task["receptor_file"]), "-l", _relativize(ligand_file)]
        if gnina_cfg.get("model_name") not in {None, "default"}:
            gnina_args.extend(["--cnn", str(gnina_cfg["model_name"])])
        command = build_gnina_command(gnina_cfg, gnina_args)
        start = time.time()
        stdout = ""
        stderr = ""
        exit_code = 127
        status = "failed"
        error = engine.warning or ""
        if engine.available:
            try:
                process = subprocess.run(
                    command,
                    check=False,
                    capture_output=True,
                    text=True,
                    timeout=int(gnina_cfg.get("timeout_seconds", 900)),
                    cwd=project_root(),
                )
                stdout = process.stdout
                stderr = process.stderr
                exit_code = process.returncode
                status = "success" if process.returncode == 0 else "failed"
                error = "" if status == "success" else stderr[-1000:]
            except Exception as exc:
                error = str(exc)
        stdout_log.write_text(stdout, encoding="utf-8")
        stderr_log.write_text(stderr or error, encoding="utf-8")
        parsed = parse_gnina_output(stdout)
        warnings = []
        if status != "success":
            warnings.append(error or "gnina_failed")
        for field in ["cnnscore", "cnnaffinity", "gnina_empirical_affinity"]:
            if parsed[field] is None:
                warnings.append(f"missing_{field}")
        runtime = round(time.time() - start, 3)
        raw_rows.append(
            {
                "rescoring_task_id": task_id,
                "pose_id": task["pose_id"],
                "gnina_version": engine.version,
                "gnina_model": task["model_name"],
                "gnina_mode": task["rescoring_mode"],
                "command_line": " ".join(command),
                "runtime_seconds": runtime,
                "exit_code": exit_code,
                "status": status,
                "stdout_log": str(stdout_log),
                "stderr_log": str(stderr_log),
                "error_message": error,
            }
        )
        score_rows.append(
            {
                "pose_id": task["pose_id"],
                "rescoring_task_id": task_id,
                "docking_task_id": task["docking_task_id"],
                "ligand_id": task["ligand_id"],
                "target_receptor_id": task["target_receptor_id"],
                "task_type": task["task_type"],
                "docking_engine": task["docking_engine"],
                "original_pose_rank": task["original_pose_rank"],
                "original_docking_score": task["original_docking_score"],
                "rmsd_symmetry_corrected": task["rmsd_symmetry_corrected"],
                "sanity_status": task["sanity_status"],
                "stage3_pose_label": task["native_like_label_stage3"],
                "gnina_version": engine.version,
                "gnina_model": task["model_name"],
                "gnina_mode": task["rescoring_mode"],
                "gnina_empirical_affinity": parsed["gnina_empirical_affinity"],
                "cnnscore": parsed["cnnscore"],
                "cnnaffinity": parsed["cnnaffinity"],
                "cnn_vs": parsed["cnn_vs"],
                "cnnscore_rank_within_task": None,
                "cnnaffinity_rank_within_task": None,
                "gnina_affinity_rank_within_task": None,
                "rescoring_status": status,
                "rescoring_warnings_json": json.dumps(warnings),
            }
        )
        raw = pd.DataFrame(raw_rows, columns=RAW_GNINA_COLUMNS).drop_duplicates("rescoring_task_id", keep="last")
        scores = pd.DataFrame(score_rows, columns=GNINA_SCORE_COLUMNS).drop_duplicates("rescoring_task_id", keep="last")
        _write_gnina_outputs(raw, scores, paths)
    raw = pd.DataFrame(raw_rows, columns=RAW_GNINA_COLUMNS).drop_duplicates("rescoring_task_id", keep="last")
    scores = pd.DataFrame(score_rows, columns=GNINA_SCORE_COLUMNS).drop_duplicates("rescoring_task_id", keep="last")
    _write_gnina_outputs(raw, scores, paths)
    return raw, scores
