"""GPU-batched Uni-Dock docking for the enrichment campaign.

Uni-Dock's ``--gpu_batch`` docks a whole list of ligands against one receptor in a
single GPU call, auto-sizing the batch to available VRAM (``--max_gpu_memory 0``).
This replaces the serial, one-subprocess-per-ligand, ``--cpu 4`` path in
``stage3/unidock_runner.py`` (which does not use the GPU at all): for a campaign of
~10^4 ligand-receptor tasks that is the difference between weeks and hours.

Design:
  * one Uni-Dock ``--gpu_batch`` invocation per receptor over ALL its ligands;
  * Uni-Dock writes ``<ligand_stem>_out.pdbqt`` per ligand into an output dir;
  * we parse each output's top pose score (best VINA RESULT) and pose file;
  * idempotent + resumable: a receptor whose per-ligand outputs already exist is
    skipped, and only missing ligands are (re)docked, so an interrupted campaign
    resumes cheaply.

This module is intentionally separate from the validated serial stage3/stage8
runners so existing results are not disturbed.
"""
from __future__ import annotations

import subprocess
import time
from pathlib import Path
from typing import Any

import pandas as pd


def _best_score_from_pdbqt(path: Path) -> float | None:
    """Return the best (most negative) Vina score across MODELs in a Uni-Dock
    output PDBQT. Uni-Dock writes 'REMARK VINA RESULT: <score> ...' per model,
    ranked best-first, so the first is the top pose."""
    if not path.exists():
        return None
    try:
        with path.open("r", encoding="utf-8", errors="replace") as fh:
            for line in fh:
                if line.startswith("REMARK VINA RESULT:"):
                    parts = line.split()
                    try:
                        return float(parts[3])
                    except (IndexError, ValueError):
                        return None
    except OSError:
        return None
    return None


def _split_top_pose(out_pdbqt: Path, dest: Path) -> bool:
    """Write only the first MODEL (top-ranked pose) of a Uni-Dock output to dest.
    Returns True on success. Used to hand a single best pose to GNINA/ProLIF."""
    if not out_pdbqt.exists():
        return False
    lines: list[str] = []
    started = False
    try:
        with out_pdbqt.open("r", encoding="utf-8", errors="replace") as fh:
            for line in fh:
                if line.startswith("MODEL"):
                    if started:  # reached the 2nd model -> stop
                        break
                    started = True
                    lines = [line]
                elif line.startswith("ENDMDL"):
                    lines.append(line)
                    break
                elif started:
                    lines.append(line)
    except OSError:
        return False
    if not lines:
        # single-model file without MODEL records: copy whole file
        try:
            dest.write_text(out_pdbqt.read_text(encoding="utf-8", errors="replace"), encoding="utf-8")
            return True
        except OSError:
            return False
    dest.write_text("".join(lines), encoding="utf-8")
    return True


def dock_receptor_batch(
    unidock: str,
    receptor_pdbqt: str | Path,
    ligand_pdbqts: list[Path],
    box: dict[str, float],
    out_dir: Path,
    *,
    search_mode: str = "balance",
    num_modes: int = 9,
    seed: int = 807,
    max_gpu_memory: int = 0,
    chunk_size: int | None = None,
    timeout_s: int = 7200,
    env_path_extra: str = "/usr/lib/wsl/lib",
) -> dict[str, Any]:
    """Dock a list of prepared PDBQT ligands against ONE receptor via a single (or
    few) Uni-Dock ``--gpu_batch`` GPU call(s). Returns a summary dict; per-ligand
    outputs land in ``out_dir`` as ``<ligand_stem>_out.pdbqt``.

    Resumable: ligands whose ``<stem>_out.pdbqt`` already exists are not re-docked.

    Ligand paths are passed to Uni-Dock via ``--ligand_index <file>`` (a text file
    of one path per line), NOT as ``--gpu_batch <path> <path> ...`` on the command
    line: a full DUD-E set is ~35k ligands and putting every path on argv exceeds
    the OS ``ARG_MAX`` (~2 MB) and hangs/crashes the exec. We still chunk (default
    ~4000/call) so a single Uni-Dock invocation's peak host memory for parsing is
    bounded and each chunk is independently resumable; Uni-Dock bins each chunk to
    VRAM internally via ``--max_gpu_memory 0``.
    """
    import os

    out_dir.mkdir(parents=True, exist_ok=True)
    env = os.environ.copy()
    env["PATH"] = env.get("PATH", "") + ":" + env_path_extra

    todo = [lig for lig in ligand_pdbqts if not (out_dir / f"{Path(lig).stem}_out.pdbqt").exists()]
    already = len(ligand_pdbqts) - len(todo)
    if not todo:
        return {"receptor": str(receptor_pdbqt), "docked": 0, "reused": already, "status": "all_cached", "calls": 0}

    cs = chunk_size if (chunk_size and chunk_size > 0) else 4000
    chunks = [todo[i : i + cs] for i in range(0, len(todo), cs)]

    calls = 0
    failures: list[str] = []
    t0 = time.time()
    for ci, chunk in enumerate(chunks):
        # write this chunk's ligand paths to an index file (avoids ARG_MAX)
        index_file = out_dir / f"_ligand_index_{ci:04d}.txt"
        index_file.write_text("\n".join(str(p) for p in chunk) + "\n", encoding="utf-8")
        cmd = [
            str(unidock),
            "--receptor", str(receptor_pdbqt),
            "--ligand_index", str(index_file),
            "--center_x", str(box["cx"]), "--center_y", str(box["cy"]), "--center_z", str(box["cz"]),
            "--size_x", str(box["sx"]), "--size_y", str(box["sy"]), "--size_z", str(box["sz"]),
            "--search_mode", search_mode,
            "--num_modes", str(num_modes),
            "--seed", str(seed),
            "--max_gpu_memory", str(max_gpu_memory),
            "--dir", str(out_dir),
        ]
        try:
            proc = subprocess.run(cmd, check=False, capture_output=True, text=True, timeout=timeout_s, env=env)
            calls += 1
            if proc.returncode != 0:
                failures.append(f"chunk{ci}: " + (proc.stderr.strip()[-400:] or proc.stdout.strip()[-400:]))
        except subprocess.TimeoutExpired:
            failures.append(f"chunk{ci}: timeout after {timeout_s}s on {len(chunk)} ligands")
        except Exception as exc:  # noqa: BLE001
            failures.append(f"chunk{ci}: {exc}")

    docked = sum(1 for lig in todo if (out_dir / f"{Path(lig).stem}_out.pdbqt").exists())
    return {
        "receptor": str(receptor_pdbqt),
        "docked": docked,
        "reused": already,
        "requested": len(todo),
        "calls": calls,
        "runtime_s": round(time.time() - t0, 1),
        "status": "ok" if docked == len(todo) else "partial",
        "failures": failures[:5],
    }


def collect_docking_scores(ligand_pdbqts: list[Path], out_dir: Path, receptor_id: str) -> pd.DataFrame:
    """Gather top-pose docking scores for a receptor's docked ligands."""
    rows = []
    for lig in ligand_pdbqts:
        stem = Path(lig).stem
        out_pdbqt = out_dir / f"{stem}_out.pdbqt"
        rows.append({
            "ligand_stem": stem,
            "target_receptor_id": receptor_id,
            "docking_score": _best_score_from_pdbqt(out_pdbqt),
            "pose_file": str(out_pdbqt) if out_pdbqt.exists() else None,
            "docked_flag": out_pdbqt.exists(),
        })
    return pd.DataFrame(rows)
