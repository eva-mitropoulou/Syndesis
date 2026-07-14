"""Batched GNINA CNN rescoring for the enrichment campaign.

The Stage-8 rescorer launches one ``docker run`` per pose; each container launch
costs ~1-2 s, so for 10^4 poses that overhead alone dominates. GNINA can score many
ligands in a single invocation when handed a multi-molecule SDF (``-l poses.sdf``),
so we group a receptor's top poses into one SDF and score them in ONE container
launch on the GPU. That collapses thousands of container starts into ~one per
receptor.

We reuse the exact score regexes from stage4/score_parser, but must split GNINA's
concatenated per-molecule output back to individual ligands (GNINA prints a score
block per input molecule, in input order).
"""
from __future__ import annotations

import re
import subprocess
from pathlib import Path
from typing import Any

import pandas as pd

from syndesis.stage4.score_parser import parse_gnina_output


def _pose_to_sdf_block(args: tuple[str, str, str, str]) -> tuple[str, str | None]:
    """Worker: convert one pose PDBQT -> SDF text block, titled with the ligand id.
    Returns (ligand_id, block_or_None). Module-level for ProcessPool pickling."""
    ligand_id, pdbqt, obabel, env_path_extra = args
    import os
    from pathlib import Path as _P

    if not _P(pdbqt).exists():
        return (ligand_id, None)
    env = os.environ.copy()
    env["PATH"] = env.get("PATH", "") + ":" + env_path_extra
    tmp_sdf = _P(pdbqt).with_suffix(".gnina_tmp.sdf")
    try:
        proc = subprocess.run(
            [obabel, str(pdbqt), "-osdf", "-O", str(tmp_sdf)],
            check=False, capture_output=True, text=True, timeout=60, env=env,
        )
    except Exception:  # noqa: BLE001
        return (ligand_id, None)
    if proc.returncode != 0 or not tmp_sdf.exists() or tmp_sdf.stat().st_size == 0:
        tmp_sdf.unlink(missing_ok=True)
        return (ligand_id, None)
    text = tmp_sdf.read_text(encoding="utf-8", errors="replace")
    tmp_sdf.unlink(missing_ok=True)
    lines = text.splitlines()
    if lines:
        lines[0] = ligand_id  # title = ligand id for re-association
    return (ligand_id, "\n".join(lines))


def build_multi_sdf(pose_pdbqts: list[tuple[str, Path]], obabel: str, out_sdf: Path,
                    env_path_extra: str = "/usr/lib/wsl/lib", max_workers: int = 10) -> list[str]:
    """Convert each top-pose PDBQT to SDF (bond orders for GNINA) IN PARALLEL and
    concatenate into one multi-molecule SDF, preserving input order (== GNINA output
    order). Serial obabel over ~35k poses took days; a process pool cuts it to
    minutes. Ligands that fail conversion are skipped and omitted from the returned
    order list."""
    from concurrent.futures import ProcessPoolExecutor

    out_sdf.parent.mkdir(parents=True, exist_ok=True)
    jobs = [(lid, str(p), obabel, env_path_extra) for lid, p in pose_pdbqts]
    results: dict[str, str] = {}
    with ProcessPoolExecutor(max_workers=max_workers) as pool:
        for lid, block in pool.map(_pose_to_sdf_block, jobs, chunksize=32):
            if block is not None:
                results[lid] = block
    # preserve original input order for GNINA output alignment
    order = [lid for lid, _ in pose_pdbqts if lid in results]
    if order:
        out_sdf.write_text("\n".join(results[lid] for lid in order) + "\n", encoding="utf-8")
    return order


def _split_gnina_per_molecule(stdout: str, n_expected: int) -> list[dict[str, float | None]]:
    """Split GNINA stdout into per-molecule score dicts.

    Verified GNINA v1.3.3 --score_only format prints, per input molecule and in
    input order, a block:
        Affinity: -9.34250 (kcal/mol)
        CNNscore: 0.11189
        CNNaffinity: 4.44162
    We split at each 'Affinity:' line and parse each block with the shared parser.
    If the block count does not match the expected molecule count, we return the
    parsed blocks we have and pad with None (so a partial parse never silently
    misaligns ligand ids to scores)."""
    blocks = re.split(r"(?im)(?=^\s*Affinity:\s)", stdout)
    blocks = [b for b in blocks if re.search(r"(?im)^\s*Affinity:\s", b)]
    parsed = [parse_gnina_output(b) for b in blocks]
    if len(parsed) < n_expected:
        parsed = parsed + [{"gnina_empirical_affinity": None, "cnnscore": None, "cnnaffinity": None, "cnn_vs": None}
                           for _ in range(n_expected - len(parsed))]
    return parsed[:n_expected]


def rescore_receptor_batch(
    docker: str,
    gnina_image: str,
    receptor_pdbqt: str | Path,
    pose_pdbqts: list[tuple[str, Path]],
    obabel: str,
    work_root: Path,
    receptor_id: str,
    *,
    cnn_scoring: str = "rescore",
    timeout_s: int = 3600,
) -> pd.DataFrame:
    """Rescore a receptor's top poses with GNINA in a single container launch.

    ``pose_pdbqts`` = [(ligand_id, top_pose_pdbqt_path), ...]. The receptor and the
    combined ligand SDF are bind-mounted into the container. Returns per-ligand
    GNINA scores. work_root should be a directory that is (or is under) a path the
    container can mount; we mount work_root itself at /work.
    """
    work_root = Path(work_root)
    work_root.mkdir(parents=True, exist_ok=True)
    multi_sdf = work_root / f"{receptor_id}_poses.sdf"
    order = build_multi_sdf(pose_pdbqts, obabel, multi_sdf)
    if not order:
        return pd.DataFrame(columns=["lit_pcba_id", "target_receptor_id", "cnnscore", "cnnaffinity", "gnina_empirical_affinity", "rescore_status"])

    # Receptor must also be reachable inside the container: copy it under work_root.
    rec_local = work_root / Path(receptor_pdbqt).name
    if not rec_local.exists():
        rec_local.write_text(Path(receptor_pdbqt).read_text(encoding="utf-8", errors="replace"), encoding="utf-8")

    cmd = [
        docker, "run", "--rm", "--gpus", "all",
        "-v", f"{work_root}:/work", "-w", "/work",
        gnina_image, "gnina",
        "--score_only",
        "-r", f"/work/{rec_local.name}",
        "-l", f"/work/{multi_sdf.name}",
        "--cnn_scoring", cnn_scoring,
    ]
    try:
        proc = subprocess.run(cmd, check=False, capture_output=True, text=True, timeout=timeout_s)
        stdout = proc.stdout
        status_ok = proc.returncode == 0
        err = proc.stderr.strip()[-500:]
    except subprocess.TimeoutExpired:
        stdout, status_ok, err = "", False, f"gnina timeout {timeout_s}s"
    except Exception as exc:  # noqa: BLE001
        stdout, status_ok, err = "", False, str(exc)

    parsed = _split_gnina_per_molecule(stdout, len(order)) if status_ok else \
        [{"gnina_empirical_affinity": None, "cnnscore": None, "cnnaffinity": None, "cnn_vs": None} for _ in order]

    rows = []
    for ligand_id, scores in zip(order, parsed, strict=False):
        rows.append({
            "lit_pcba_id": ligand_id,
            "target_receptor_id": receptor_id,
            "cnnscore": scores.get("cnnscore"),
            "cnnaffinity": scores.get("cnnaffinity"),
            "gnina_empirical_affinity": scores.get("gnina_empirical_affinity"),
            "rescore_status": "success" if status_ok else f"failed:{err}",
        })
    return pd.DataFrame(rows)
