"""Enrichment campaign orchestrator: LIT-PCBA actives+decoys across the ensemble.

For each ensemble receptor, run the GPU-saturating pipeline:
  1. Uni-Dock ``--gpu_batch`` over ALL ligands for this receptor (one GPU call);
  2. split each ligand's top pose;
  3. GNINA batch CNN rescore (one container for the receptor);
  4. ProLIF interaction recall vs the native consensus fingerprint (parallel, CPU),
     which supplies the "+interaction-constraint" signal.

The result is one master table (ligand x receptor x scores) that feeds three matched
scoring arms downstream: docking-score-only, +interaction-constraint, +pose-confidence.

Resumable: each receptor's outputs are checkpointed to parquet; a completed receptor
is skipped on re-run. Docking is itself resumable per-ligand (see gpu_docking).
Receptors are processed sequentially (the GPU is the bottleneck) but each receptor's
docking (GPU) and the previous receptor's ProLIF (CPU) can overlap via a worker.
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pandas as pd

from egfr_dockingforge.enrichment.gpu_docking import (
    collect_docking_scores,
    dock_receptor_batch,
    _split_top_pose,
)
from egfr_dockingforge.enrichment.gnina_batch import rescore_receptor_batch


def _box_for_receptor(receptors: pd.DataFrame, receptor_id: str) -> dict[str, float]:
    import ast

    row = receptors[receptors["receptor_id"] == receptor_id].iloc[0]

    def _vec(v):
        return [float(x) for x in (ast.literal_eval(v) if isinstance(v, str) else v)]

    c = _vec(row["suggested_docking_box_center"])
    s = _vec(row["suggested_docking_box_size"])
    return {"cx": c[0], "cy": c[1], "cz": c[2], "sx": s[0], "sy": s[1], "sz": s[2]}


def run_receptor(
    receptor_id: str,
    receptor_pdbqt: str,
    box: dict[str, float],
    ligands: pd.DataFrame,   # columns: lit_pcba_id, prepared_pdbqt, label, activity
    tools: dict[str, str],
    work_root: Path,
    *,
    num_modes: int = 9,
    seed: int = 807,
) -> pd.DataFrame:
    """Dock + rescore all ligands against one receptor. Returns the per-ligand
    score table for this receptor (checkpointed by the caller)."""
    rec_work = work_root / receptor_id
    dock_dir = rec_work / "dock"
    top_dir = rec_work / "top"
    top_dir.mkdir(parents=True, exist_ok=True)

    lig_paths = [Path(p) for p in ligands["prepared_pdbqt"].tolist() if p and Path(p).exists()]

    # 1. GPU batch dock (resumable inside).
    dock_summary = dock_receptor_batch(
        tools["unidock"], receptor_pdbqt, lig_paths, box, dock_dir,
        num_modes=num_modes, seed=seed,
    )

    # 2. Collect docking scores + split top poses for the poses that docked.
    scores = collect_docking_scores(lig_paths, dock_dir, receptor_id)
    id_by_stem = {Path(p).stem: lid for lid, p in zip(ligands["lit_pcba_id"], ligands["prepared_pdbqt"]) if p}
    scores["lit_pcba_id"] = scores["ligand_stem"].map(id_by_stem)

    pose_pairs: list[tuple[str, Path]] = []
    for row in scores[scores["docked_flag"]].to_dict("records"):
        top = top_dir / f"{row['ligand_stem']}_top.pdbqt"
        if not top.exists():
            _split_top_pose(Path(row["pose_file"]), top)
        if top.exists():
            pose_pairs.append((row["lit_pcba_id"], top))

    # 3. GNINA batch rescore (one container for this receptor).
    gnina = rescore_receptor_batch(
        tools["docker"], tools["gnina_image"], receptor_pdbqt, pose_pairs,
        tools["obabel"], rec_work / "gnina", receptor_id,
    )

    # 4. Merge docking + GNINA into the receptor score table. (ProLIF interaction
    #    recall is added in a later step that reuses the stage8 ProLIF engine on the
    #    top-pose SDFs; kept separate so a GNINA-only pass is still useful.)
    merged = scores.merge(ligands[["lit_pcba_id", "label", "activity"]], on="lit_pcba_id", how="left")
    if not gnina.empty:
        merged = merged.merge(
            gnina[["lit_pcba_id", "cnnscore", "cnnaffinity", "gnina_empirical_affinity", "rescore_status"]],
            on="lit_pcba_id", how="left",
        )
    else:
        for c in ["cnnscore", "cnnaffinity", "gnina_empirical_affinity"]:
            merged[c] = None
        merged["rescore_status"] = "no_poses"
    merged["target_receptor_id"] = receptor_id
    merged["dock_summary_json"] = json.dumps(dock_summary)
    return merged


def run_campaign(
    ligands: pd.DataFrame,
    receptors: pd.DataFrame,
    tools: dict[str, str],
    work_root: Path,
    *,
    receptor_ids: list[str] | None = None,
    num_modes: int = 9,
    seed: int = 807,
) -> pd.DataFrame:
    """Run the full ligands x receptors campaign, checkpointing per receptor.
    Returns the concatenated master score table across receptors."""
    work_root = Path(work_root)
    work_root.mkdir(parents=True, exist_ok=True)
    if receptor_ids is None:
        receptor_ids = receptors[receptors["selected_flag"].fillna(False).astype(bool)]["receptor_id"].tolist()

    rec_prep = receptors.set_index("receptor_id")
    all_tables = []
    for rid in receptor_ids:
        ckpt = work_root / f"scores_{rid}.parquet"
        if ckpt.exists():
            all_tables.append(pd.read_parquet(ckpt))
            continue
        receptor_pdbqt = _resolve_receptor_pdbqt(rec_prep.loc[rid], tools)
        box = _box_for_receptor(receptors, rid)
        table = run_receptor(rid, receptor_pdbqt, box, ligands, tools, work_root,
                             num_modes=num_modes, seed=seed)
        table.to_parquet(ckpt, index=False)
        all_tables.append(table)

    master = pd.concat(all_tables, ignore_index=True) if all_tables else pd.DataFrame()
    master.to_parquet(work_root / "enrichment_master_scores.parquet", index=False)
    master.to_csv(work_root / "enrichment_master_scores.csv", index=False)
    return master


def _resolve_receptor_pdbqt(rec_row: pd.Series, tools: dict[str, str]) -> str:
    """Prefer the stage3-prepared docking pdbqt for the receptor; fall back to a
    configured docking_receptors dir keyed by receptor_id."""
    for key in ("docking_format_file", "prepared_receptor_file"):
        val = rec_row.get(key) if hasattr(rec_row, "get") else None
        if val and str(val).endswith(".pdbqt") and Path(val).exists():
            return str(val)
    # fall back: <docking_receptors>/<receptor_id>.pdbqt
    cand = Path(tools["docking_receptors_dir"]) / f"{rec_row.name}.pdbqt"
    if cand.exists():
        return str(cand)
    raise FileNotFoundError(f"No prepared receptor pdbqt for {rec_row.name}")
