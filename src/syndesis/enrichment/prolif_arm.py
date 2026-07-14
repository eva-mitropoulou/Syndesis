"""ProLIF interaction-recall for enrichment poses (the +interaction-constraint arm).

For each docked top pose, compute its ProLIF interaction fingerprint against the
docked receptor and measure recall/Tanimoto vs the native consensus fingerprint
(the union of key interactions the crystallographic ligands make). This is the
signal that distinguishes a pose that genuinely recovers the validated EGFR binding
mode from a docking-score "hacker" that scores well without the right contacts.

Reuses the validated stage5 ProLIF engine. Embarrassingly parallel (CPU-bound), so
it runs across a process pool on the saved pose files after GPU docking finishes.
Idempotent: writes a per-receptor parquet checkpoint and skips completed receptors.
"""
from __future__ import annotations

import json
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Any

import pandas as pd


def _consensus_bits(native_fps: pd.DataFrame) -> set[str]:
    """Union of native interaction-fingerprint bits across the crystallographic
    complexes = the consensus binding-mode target."""
    bits: set[str] = set()
    for v in native_fps["fingerprint_sparse_json"].dropna():
        payload = json.loads(v) if isinstance(v, str) else v
        bits.update(payload)
    return bits


def _posed_ligand_sdf(pose_pdbqt: str, template_sdf: str, lig_dir: str, _obabel: str) -> str:
    """Reconstruct the docked pose using the prepared, bond-order-valid graph."""
    from pathlib import Path as _P
    from syndesis.stage5.pose_reconstruction import reconstruct_pose_sdf

    out_sdf = _P(lig_dir) / f"{_P(pose_pdbqt).stem}.pose_template.sdf"
    _P(lig_dir).mkdir(parents=True, exist_ok=True)
    if out_sdf.exists() and out_sdf.stat().st_size > 0:
        return str(out_sdf)
    prepared_pdbqt = str(_P(template_sdf).with_suffix(".pdbqt"))
    return str(reconstruct_pose_sdf(pose_pdbqt, template_sdf, prepared_pdbqt, out_sdf))


def _score_one_pose(args: tuple) -> dict[str, Any]:
    """Worker: compute interaction recall for one pose. Self-contained for pickling.
    args = (ligand_id, receptor_id, receptor_prolif_pdb, pose_pdbqt, prolif_lig_dir,
            residue_map_records, config, consensus_bits_list, template_sdf, obabel)"""
    (ligand_id, receptor_id, receptor_pdb, pose_pdbqt, lig_dir,
     residue_map_records, config, consensus_list, template_sdf, obabel) = args
    try:
        from syndesis.stage5.prolif_engine import fingerprint_from_interactions, tanimoto
        from syndesis.stage8.screening_interactions import _compute_interactions_rdkit
        import pandas as _pd

        consensus = set(consensus_list)
        residue_map = _pd.DataFrame(residue_map_records)
        lig_sdf = _posed_ligand_sdf(pose_pdbqt, template_sdf, lig_dir, obabel)
        interactions, _meta = _compute_interactions_rdkit(
            Path(receptor_pdb), Path(lig_sdf), residue_map, config
        )
        _bs, _js, bits = fingerprint_from_interactions(interactions)
        recovered = bits & consensus
        recall = len(recovered) / len(consensus) if consensus else 0.0
        return {
            "lit_pcba_id": ligand_id,
            "target_receptor_id": receptor_id,
            "num_interactions": len(bits),
            "key_interaction_recall_consensus": recall,
            "ifp_tanimoto_to_consensus": tanimoto(bits, consensus) or 0.0,
            "prolif_status": "ok",
        }
    except Exception as exc:  # noqa: BLE001
        return {
            "lit_pcba_id": ligand_id,
            "target_receptor_id": receptor_id,
            "num_interactions": 0,
            "key_interaction_recall_consensus": 0.0,
            "ifp_tanimoto_to_consensus": 0.0,
            "prolif_status": f"failed:{type(exc).__name__}:{str(exc)[:120]}",
        }


def score_receptor_interactions(
    receptor_id: str,
    receptor_prolif_pdb: Path,
    pose_pairs: list[tuple[str, Path]],   # (ligand_id, top_pose_pdbqt)
    native_fps: pd.DataFrame,
    residue_map: pd.DataFrame,
    config: dict[str, Any],
    work_root: Path,
    prepared_dir: Path,
    obabel: str,
    *,
    max_workers: int = 10,
) -> pd.DataFrame:
    """Compute interaction recall for all of a receptor's docked top poses, in
    parallel. Checkpointed to interactions_<receptor>.parquet (resumable).
    prepared_dir holds the per-ligand template SDFs (<lit_pcba_id>.sdf) used to
    rebuild chemically-correct, H-complete ligands at the docked coordinates."""
    ckpt = work_root / f"interactions_{receptor_id}.parquet"
    if ckpt.exists():
        return pd.read_parquet(ckpt)
    lig_dir = work_root / receptor_id / "prolif_ligands"
    lig_dir.mkdir(parents=True, exist_ok=True)
    consensus = list(_consensus_bits(native_fps))
    residue_records = residue_map.to_dict("records")
    jobs = [
        (lid, receptor_id, str(receptor_prolif_pdb), str(pose), str(lig_dir),
         residue_records, config, consensus, str(Path(prepared_dir) / f"{lid}.sdf"), obabel)
        for lid, pose in pose_pairs
    ]
    rows: list[dict[str, Any]] = []
    with ProcessPoolExecutor(max_workers=max_workers) as pool:
        futs = [pool.submit(_score_one_pose, j) for j in jobs]
        for fut in as_completed(futs):
            rows.append(fut.result())
    out = pd.DataFrame(rows)
    out.to_parquet(ckpt, index=False)
    failures = out.loc[~out["prolif_status"].eq("ok"), ["lit_pcba_id", "prolif_status"]]
    if not failures.empty:
        raise RuntimeError(
            f"{receptor_id}: {len(failures)} ProLIF calculations failed; "
            f"ranking is blocked. Examples: {failures.head(5).to_dict('records')}"
        )
    return out


def add_interaction_arm(
    master: pd.DataFrame,
    campaign_root: Path,
    receptors: pd.DataFrame,
    native_fps: pd.DataFrame,
    residue_map: pd.DataFrame,
    config: dict[str, Any],
    work_root: Path,
    prepared_dir: Path,
    obabel: str,
    *,
    max_workers: int = 10,
) -> pd.DataFrame:
    """For each receptor in the master table, locate its docked top poses and
    compute ProLIF interaction recall, then merge onto master so build_arms() can
    form the real arm_inter. Reuses top-pose PDBQTs produced by the orchestrator."""
    from syndesis.stage5.prolif_engine import prepare_protein_for_prolif

    rec_prep = receptors.set_index("receptor_id")
    all_int = []
    for rid, grp in master.groupby("target_receptor_id"):
        top_dir = campaign_root / rid / "top"
        pose_pairs = []
        for row in grp.to_dict("records"):
            stem = row.get("ligand_stem")
            if not stem:
                continue
            top = top_dir / f"{stem}_top.pdbqt"
            if top.exists():
                pose_pairs.append((row["lit_pcba_id"], top))
        if not pose_pairs:
            continue
        # receptor pdb for ProLIF (explicit-H); reuse stage3 prepared receptor pdb
        rec_row = rec_prep.loc[rid]
        rec_src = rec_row.get("receptor_file_path") or rec_row.get("aligned_receptor_file_path")
        rec_pdb = prepare_protein_for_prolif(rec_src, work_root / rid / "prolif_receptor")
        idf = score_receptor_interactions(rid, rec_pdb, pose_pairs, native_fps, residue_map, config, work_root, prepared_dir, obabel, max_workers=max_workers)
        all_int.append(idf)
    if not all_int:
        return master
    inter = pd.concat(all_int, ignore_index=True)
    merged = master.merge(
        inter[["lit_pcba_id", "target_receptor_id", "key_interaction_recall_consensus", "ifp_tanimoto_to_consensus"]],
        on=["lit_pcba_id", "target_receptor_id"], how="left",
    )
    merged.to_parquet(work_root / "enrichment_master_scores_with_interactions.parquet", index=False)
    return merged
