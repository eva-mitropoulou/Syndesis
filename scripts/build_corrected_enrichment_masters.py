"""Build fail-closed enrichment masters from strict pose-graph ProLIF bitsets."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
ROBUSTNESS = ROOT / "data/processed/submission_robustness"
OUTPUT = ROOT / "results/analysis_inputs"
EXTERNAL = ROOT / "data/external"

TARGETS = {
    "egfr": {
        "master": EXTERNAL / "egfr/enrichment_master_scores_with_interactions.parquet",
        "native": EXTERNAL / "egfr/native_interaction_fingerprints.parquet",
    },
    "cdk2": {
        "master": EXTERNAL / "cdk2/enrichment_master_scores_with_interactions.parquet",
        "native": EXTERNAL / "cdk2/native_interaction_fingerprints.parquet",
    },
}


def parse_bits(value: str) -> set[str]:
    return set(json.loads(value))


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--target", choices=sorted(TARGETS), action="append")
    args = parser.parse_args()
    OUTPUT.mkdir(parents=True, exist_ok=True)
    selected = args.target or list(TARGETS)
    for target in selected:
        paths = TARGETS[target]
        fingerprint_paths = sorted((ROOT / "data/processed/pose_fingerprints" / target).glob("pose_fingerprints_*.parquet"))
        if len(fingerprint_paths) != 5:
            raise RuntimeError(f"{target}: expected five fingerprint checkpoints, found {len(fingerprint_paths)}")
        fingerprints = pd.concat((pd.read_parquet(path) for path in fingerprint_paths), ignore_index=True)
        if fingerprints.duplicated(["ligand_id", "target_receptor_id"]).any():
            raise RuntimeError(f"{target}: duplicate strict fingerprint rows")
        fingerprints["lit_pcba_id"] = fingerprints["ligand_id"].astype(str)

        master = pd.read_parquet(paths["master"]).sort_values("cnnscore", ascending=False).drop_duplicates(
            ["lit_pcba_id", "target_receptor_id"]
        )
        master["lit_pcba_id"] = master["lit_pcba_id"].astype(str)
        drop = ["key_interaction_recall_consensus", "ifp_tanimoto_to_consensus"]
        master = master.drop(columns=[column for column in drop if column in master], errors="ignore")
        corrected = master.merge(
            fingerprints[[
                "lit_pcba_id", "target_receptor_id", "fingerprint_sparse_json", "num_interactions", "status", "error",
            ]],
            on=["lit_pcba_id", "target_receptor_id"],
            how="left",
            validate="one_to_one",
        )
        missing = corrected["status"].isna()
        if missing.any():
            raise RuntimeError(f"{target}: {int(missing.sum())} score rows lack strict fingerprints")

        failed = ~corrected["status"].eq("ok")
        empty_pose = corrected["error"].fillna("").str.contains("posed=0", regex=False)
        no_scored_pose = failed & empty_pose & corrected["cnnscore"].isna()
        corrected.loc[no_scored_pose, "status"] = "no_scored_pose"
        corrected.loc[no_scored_pose, "fingerprint_sparse_json"] = None
        corrected.loc[no_scored_pose, "num_interactions"] = pd.NA
        failures = corrected[~corrected["status"].isin(["ok", "no_scored_pose"])]
        if not failures.empty:
            failures.to_parquet(OUTPUT / f"{target}_fingerprint_failures.parquet", index=False)
            raise RuntimeError(
                f"{target}: {len(failures)} strict ProLIF calculations failed for scored poses; corrected ranking is blocked"
            )

        native = pd.read_parquet(paths["native"])
        native_union: set[str] = set()
        for value in native["fingerprint_sparse_json"]:
            native_union.update(parse_bits(value))
        valid = corrected["status"].eq("ok")
        corrected.loc[valid, "pose_bits"] = corrected.loc[valid, "fingerprint_sparse_json"].map(parse_bits)
        corrected.loc[valid, "key_interaction_recall_native_union"] = corrected.loc[valid, "pose_bits"].map(
            lambda value: len(value & native_union) / len(native_union)
        )
        corrected.loc[valid, "ifp_jaccard_native_union"] = corrected.loc[valid, "pose_bits"].map(
            lambda value: len(value & native_union) / len(value | native_union) if value | native_union else 0.0
        )
        corrected = corrected.drop(columns=["pose_bits"])
        corrected = corrected.rename(columns={
            "key_interaction_recall_native_union": "key_interaction_recall_consensus",
            "ifp_jaccard_native_union": "ifp_tanimoto_to_consensus",
        })
        corrected.to_parquet(OUTPUT / f"{target}_master.parquet", index=False)
        metadata = {
            "target": target.upper(),
            "n_rows": len(corrected),
            "n_ligands": corrected["lit_pcba_id"].nunique(),
            "n_receptors": corrected["target_receptor_id"].nunique(),
            "native_union_size": len(native_union),
            "fingerprint_failures": 0,
            "no_scored_pose_rows": int(no_scored_pose.sum()),
            "pose_graph_policy": "Prepared SDF molecular graph with strictly mapped docked coordinates",
        }
        (OUTPUT / f"{target}_metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")
        print(json.dumps(metadata), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
