"""Recompute and archive per-pose ProLIF bitsets for native-prior sensitivity analyses."""
from __future__ import annotations

import argparse
import json
import os
import tempfile
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import pandas as pd
import yaml

from egfr_dockingforge.stage5.prolif_engine import fingerprint_from_interactions
from egfr_dockingforge.stage5.pose_reconstruction import reconstruct_pose_sdf
from egfr_dockingforge.stage8.screening_interactions import _compute_interactions_rdkit


def reconstructed_pose_sdf(ligand_id: str, posed_pdb: str, template_root: str) -> Path:
    """Transfer docked coordinates onto the original, chemically valid ligand graph."""
    root = Path(template_root)
    template_sdf = root / f"{ligand_id}.sdf"
    prepared_pdbqt = root / f"{ligand_id}.pdbqt"
    handle = tempfile.NamedTemporaryFile(prefix=f"{ligand_id}_", suffix=".sdf", delete=False)
    handle.close()
    target = Path(handle.name)
    reconstruct_pose_sdf(posed_pdb, template_sdf, prepared_pdbqt, target)
    return target


def score_pose(job: tuple[str, str, str, str, list[dict], dict]) -> dict:
    ligand_id, receptor_id, ligand_file, template_root, residue_records, config = job
    receptor_file = config.pop("_receptor_file")
    pose_sdf = None
    try:
        pose_sdf = reconstructed_pose_sdf(ligand_id, ligand_file, template_root)
        interactions, metadata = _compute_interactions_rdkit(
            Path(receptor_file),
            pose_sdf,
            pd.DataFrame(residue_records),
            config,
        )
        _bitstring, sparse_json, bits = fingerprint_from_interactions(interactions)
        return {
            "ligand_id": ligand_id,
            "target_receptor_id": receptor_id,
            "fingerprint_sparse_json": sparse_json,
            "num_interactions": len(bits),
            "interaction_engine": metadata.get("interaction_engine"),
            "interaction_engine_version": metadata.get("interaction_engine_version"),
            "status": "ok",
            "error": "",
        }
    except Exception as exc:  # noqa: BLE001
        return {
            "ligand_id": ligand_id,
            "target_receptor_id": receptor_id,
            "fingerprint_sparse_json": "[]",
            "num_interactions": 0,
            "interaction_engine": "prolif",
            "interaction_engine_version": "2.2.0",
            "status": "failed",
            "error": f"{type(exc).__name__}: {str(exc)[:240]}",
        }
    finally:
        if pose_sdf is not None:
            pose_sdf.unlink(missing_ok=True)


def ligand_id(path: Path) -> str:
    suffix = "_top.posed_h.pdb"
    if not path.name.endswith(suffix):
        raise ValueError(f"Unexpected pose filename: {path}")
    return path.name[: -len(suffix)]


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-root", type=Path, required=True)
    parser.add_argument("--output-root", type=Path, required=True)
    parser.add_argument("--residue-map", type=Path, required=True)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--template-root", type=Path, required=True)
    parser.add_argument("--workers", type=int, default=max(1, min(56, (os.cpu_count() or 2) - 4)))
    parser.add_argument("--batch-size", type=int, default=7000)
    args = parser.parse_args()

    args.output_root.mkdir(parents=True, exist_ok=True)
    residue_map = pd.read_parquet(args.residue_map)
    base_config = yaml.safe_load(args.config.read_text())
    receptor_dirs = sorted(p for p in args.input_root.iterdir() if p.is_dir())

    for receptor_dir in receptor_dirs:
        receptor_id = receptor_dir.name
        checkpoint = args.output_root / f"pose_fingerprints_{receptor_id}.parquet"
        partial = args.output_root / f"pose_fingerprints_{receptor_id}.partial.parquet"
        if checkpoint.exists():
            print(f"[{receptor_id}] checkpoint exists: {checkpoint}", flush=True)
            continue
        receptor_files = list((receptor_dir / "prolif_receptor").glob("*.h.pdb"))
        if len(receptor_files) != 1:
            raise RuntimeError(f"Expected one receptor PDB for {receptor_id}, found {len(receptor_files)}")
        poses = sorted((receptor_dir / "prolif_ligands").glob("*_top.posed_h.pdb"))
        if not poses:
            raise RuntimeError(f"No posed ligand files found for {receptor_id}")

        receptor_map = residue_map[residue_map["receptor_id"] == receptor_id]
        residue_records = receptor_map.to_dict("records")
        config = dict(base_config)
        config["_receptor_file"] = str(receptor_files[0])
        all_jobs = [
            (ligand_id(pose), receptor_id, str(pose), str(args.template_root), residue_records, dict(config))
            for pose in poses
        ]
        frame = pd.read_parquet(partial) if partial.exists() else pd.DataFrame()
        completed_ids = set(frame["ligand_id"].astype(str)) if not frame.empty else set()
        jobs = [job for job in all_jobs if job[0] not in completed_ids]
        print(
            f"[{receptor_id}] {len(all_jobs)} poses, {len(completed_ids)} resumed, "
            f"{args.workers} workers, batch size {args.batch_size}",
            flush=True,
        )
        for start in range(0, len(jobs), args.batch_size):
            batch = jobs[start : start + args.batch_size]
            rows = []
            with ProcessPoolExecutor(max_workers=args.workers) as pool:
                futures = [pool.submit(score_pose, job) for job in batch]
                for completed, future in enumerate(as_completed(futures), start=1):
                    rows.append(future.result())
                    total = len(completed_ids) + start + completed
                    if completed % 1000 == 0 or completed == len(futures):
                        print(f"[{receptor_id}] {total}/{len(all_jobs)}", flush=True)
            frame = pd.concat([frame, pd.DataFrame(rows)], ignore_index=True)
            frame = frame.drop_duplicates(["ligand_id", "target_receptor_id"]).sort_values("ligand_id")
            frame.to_parquet(partial, index=False)
            print(f"[{receptor_id}] partial checkpoint: {len(frame)}/{len(all_jobs)}", flush=True)
        if len(frame) != len(all_jobs):
            raise RuntimeError(f"{receptor_id}: expected {len(all_jobs)} rows, found {len(frame)}")
        frame.to_parquet(checkpoint, index=False)
        partial.unlink(missing_ok=True)
        status = frame["status"].value_counts(dropna=False).to_dict()
        print(f"[{receptor_id}] wrote {checkpoint}; status={json.dumps(status)}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
