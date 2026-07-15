"""Build strict ProLIF fingerprints for the corrected 1PXN docking rerun.

This consumes only the fresh 1PXN top poses produced by
``rerun_cdk2_1pxn_altloc.py``.  It reconstructs every pose on its prepared SDF
graph before interaction calculation and fails closed for every scored pose.
"""
from __future__ import annotations

import json
import os
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import pandas as pd
import yaml

from export_enrichment_pose_fingerprints import score_pose
from syndesis.stage5.prolif_engine import prepare_protein_for_prolif


ROOT = Path(__file__).resolve().parents[1]
RERUN = Path(os.environ.get("SYNTHESIS_1PXN_RERUN", "cdk2_1pxn_altloc_rerun"))
CDK2 = Path(os.environ.get("SYNTHESIS_CDK2_WORK", "cdk2_dude"))
RECEPTOR_ID = "1pxn_a_ck6"
WORKERS = max(1, min(56, (os.cpu_count() or 2) - 4))
BATCH_SIZE = 2_000


def main() -> int:
    score_path = RERUN / f"scores_{RECEPTOR_ID}.parquet"
    if not score_path.exists() or not (RERUN / "RERUN_DONE").exists():
        raise RuntimeError("The corrected 1PXN docking/GNINA rerun is not complete")

    output = RERUN / f"pose_fingerprints_{RECEPTOR_ID}.parquet"
    partial = RERUN / f"pose_fingerprints_{RECEPTOR_ID}.partial.parquet"
    scores = pd.read_parquet(score_path).copy()
    scores["lit_pcba_id"] = scores["lit_pcba_id"].astype(str)
    if len(scores) != 28_296 or scores["lit_pcba_id"].nunique() != 28_296:
        raise RuntimeError("Corrected 1PXN scores do not contain the frozen 28,296-ligand set")

    receptor_source = CDK2 / "receptors" / "1pxn_A_protein.pdb"
    receptor = prepare_protein_for_prolif(receptor_source, RERUN / "prolif_receptor")
    residue_map = pd.read_parquet(CDK2 / "stage5" / "interaction_residue_map.parquet")
    config = yaml.safe_load((ROOT / "configs/stage8_candidate_screening.yaml").read_text())
    config["_receptor_file"] = str(receptor)

    completed = pd.read_parquet(partial) if partial.exists() else pd.DataFrame()
    completed_ids = set(completed["ligand_id"].astype(str)) if not completed.empty else set()
    scored = scores[scores["cnnscore"].notna()].copy()
    rows = []
    for row in scores[scores["cnnscore"].isna()].to_dict("records"):
        rows.append({
            "ligand_id": str(row["lit_pcba_id"]),
            "target_receptor_id": RECEPTOR_ID,
            "fingerprint_sparse_json": None,
            "num_interactions": pd.NA,
            "interaction_engine": "prolif",
            "interaction_engine_version": "2.2.0",
            "status": "no_scored_pose",
            "error": "no_scored_pose",
        })

    jobs = []
    top_dir = RERUN / "campaign" / RECEPTOR_ID / "top"
    for row in scored.to_dict("records"):
        ligand = str(row["lit_pcba_id"])
        if ligand in completed_ids:
            continue
        pose = top_dir / f"{row['ligand_stem']}_top.pdbqt"
        if not pose.exists():
            raise FileNotFoundError(f"Scored 1PXN pose is missing: {pose}")
        jobs.append((ligand, RECEPTOR_ID, str(pose), str(CDK2 / "prepared_ligands"), residue_map.to_dict("records"), dict(config)))

    frame = completed.copy()
    if rows:
        frame = pd.concat([frame, pd.DataFrame(rows)], ignore_index=True)
        frame = frame.drop_duplicates(["ligand_id", "target_receptor_id"], keep="last")
        frame.to_parquet(partial, index=False)
    print(f"{len(scored)} scored poses; {len(jobs)} need strict fingerprints; {WORKERS} workers", flush=True)
    for start in range(0, len(jobs), BATCH_SIZE):
        batch = jobs[start : start + BATCH_SIZE]
        output_rows = []
        with ProcessPoolExecutor(max_workers=WORKERS) as pool:
            futures = [pool.submit(score_pose, job) for job in batch]
            for completed_count, future in enumerate(as_completed(futures), start=1):
                output_rows.append(future.result())
                if completed_count % 500 == 0 or completed_count == len(futures):
                    print(f"{start + completed_count}/{len(jobs)} newly fingerprinted", flush=True)
        frame = pd.concat([frame, pd.DataFrame(output_rows)], ignore_index=True)
        frame = frame.drop_duplicates(["ligand_id", "target_receptor_id"], keep="last")
        frame.to_parquet(partial, index=False)

    if len(frame) != len(scores):
        raise RuntimeError(f"Expected {len(scores)} fingerprint rows, found {len(frame)}")
    failures = frame[~frame["status"].isin(["ok", "no_scored_pose"])]
    if not failures.empty:
        failures.to_csv(RERUN / "strict_fingerprint_failures.csv", index=False)
        raise RuntimeError(f"{len(failures)} scored 1PXN poses failed strict fingerprinting")
    frame = frame.sort_values("ligand_id").reset_index(drop=True)
    frame.to_parquet(output, index=False)
    frame.to_csv(output.with_suffix(".csv"), index=False)
    partial.unlink(missing_ok=True)
    (RERUN / "FINGERPRINTS_DONE").write_text("done\n")
    print(json.dumps(frame["status"].value_counts().to_dict()), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
