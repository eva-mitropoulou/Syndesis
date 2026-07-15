"""Rebuild the four-receptor/four-prior CDK2 results after the 1PXN fix.

Only the corrected 1PXN docking/GNINA/strict-fingerprint rows replace their
historical counterparts.  The script writes a self-contained staging package and
does not overwrite the paper-facing result directory.
"""
from __future__ import annotations

import json
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd

import submission_robustness_analysis as statistics
from native_prior_sensitivity_analysis import bits, evaluate, jaccard, recall, weighted_recall


ROOT = Path(__file__).resolve().parents[1]
RERUN = Path("/home/dimit/cdk2_1pxn_altloc_rerun")
ROBUSTNESS = ROOT / "results_showcase" / "submission_robustness"
SOURCE_MASTER = ROBUSTNESS / "corrected_enrichment" / "cdk2_master.parquet"
NATIVE = Path("/mnt/e/cdk2_dude/stage5/native_interaction_fingerprints.parquet")
OUT = RERUN / "four_receptor_four_prior_results"
PRIMARY = ("1fin_a_atp", "2a4l_a_rrc", "1aq1_a_stu", "1pxn_a_ck6")


def parse_bits(value: str | None) -> set[str]:
    return bits(value) if isinstance(value, str) and value else set()


def primary_native_union() -> tuple[pd.DataFrame, set[str], Counter[str]]:
    native = pd.read_parquet(NATIVE).copy()
    native = native[native["receptor_id"].isin(PRIMARY)].copy()
    if set(native["receptor_id"]) != set(PRIMARY):
        raise RuntimeError("The four primary CDK2 native complexes are incomplete")
    native["native_bits"] = native["fingerprint_sparse_json"].map(parse_bits)
    counts = Counter(bit for value in native["native_bits"] for bit in value)
    union = set(counts)
    if len(union) != 38:
        raise RuntimeError(f"Expected the frozen four-native CDK2 prior to contain 38 bits, found {len(union)}")
    return native, union, counts


def build_master(union: set[str]) -> pd.DataFrame:
    old = pd.read_parquet(SOURCE_MASTER).copy()
    old["lit_pcba_id"] = old["lit_pcba_id"].astype(str)
    old = old[old["target_receptor_id"].isin(PRIMARY)].copy()
    old = old[old["target_receptor_id"] != "1pxn_a_ck6"].copy()

    scores = pd.read_parquet(RERUN / "scores_1pxn_a_ck6.parquet").copy()
    fingerprints = pd.read_parquet(RERUN / "pose_fingerprints_1pxn_a_ck6.parquet").copy()
    scores["lit_pcba_id"] = scores["lit_pcba_id"].astype(str)
    fingerprints["ligand_id"] = fingerprints["ligand_id"].astype(str)
    replaced = scores.merge(
        fingerprints[["ligand_id", "fingerprint_sparse_json", "num_interactions", "status", "error"]],
        left_on="lit_pcba_id", right_on="ligand_id", how="left", validate="one_to_one",
    ).drop(columns="ligand_id")
    if replaced["status"].isna().any():
        raise RuntimeError("Corrected 1PXN scores lack strict fingerprint rows")
    failed = ~replaced["status"].isin(["ok", "no_scored_pose"])
    if failed.any():
        raise RuntimeError(f"Corrected 1PXN has {int(failed.sum())} strict fingerprint failures")

    master = pd.concat([old, replaced], ignore_index=True, sort=False)
    master = master.sort_values("cnnscore", ascending=False, na_position="last").drop_duplicates(
        ["lit_pcba_id", "target_receptor_id"], keep="first"
    )
    if len(master) != 113_184 or master["lit_pcba_id"].nunique() != 28_296:
        raise RuntimeError(f"Unexpected four-receptor master dimensions: {master.shape}")
    counts = master["target_receptor_id"].value_counts().to_dict()
    if any(counts.get(receptor) != 28_296 for receptor in PRIMARY):
        raise RuntimeError(f"Incomplete receptor coverage: {counts}")

    valid = master["status"].eq("ok")
    master["key_interaction_recall_consensus"] = np.nan
    master["ifp_tanimoto_to_consensus"] = np.nan
    master.loc[valid, "pose_bits"] = master.loc[valid, "fingerprint_sparse_json"].map(parse_bits)
    master.loc[valid, "key_interaction_recall_consensus"] = master.loc[valid, "pose_bits"].map(
        lambda observed: recall(observed, union)
    )
    master.loc[valid, "ifp_tanimoto_to_consensus"] = master.loc[valid, "pose_bits"].map(
        lambda observed: jaccard(observed, union)
    )
    return master.drop(columns="pose_bits")


def prior_sensitivity(master: pd.DataFrame, native: pd.DataFrame, union: set[str], counts: Counter[str]) -> pd.DataFrame:
    raw = master.copy()
    raw["pose_bits"] = raw["fingerprint_sparse_json"].map(parse_bits)
    native_sets = dict(zip(native["receptor_id"], native["native_bits"]))
    core = {bit for bit, count in counts.items() if count / len(native) >= 0.60}
    weights = {bit: count / len(native) for bit, count in counts.items()}
    terms = {
        "primary_union_recall": raw["pose_bits"].map(lambda observed: recall(observed, union)),
        "60pct_core_recall": raw["pose_bits"].map(lambda observed: recall(observed, core)),
        "frequency_weighted_recall": raw["pose_bits"].map(lambda observed: weighted_recall(observed, weights)),
        "primary_union_jaccard": raw["pose_bits"].map(lambda observed: jaccard(observed, union)),
        "receptor_specific_recall": pd.Series([
            recall(observed, native_sets[receptor]) if receptor in native_sets else np.nan
            for observed, receptor in zip(raw["pose_bits"], raw["target_receptor_id"])
        ], index=raw.index),
    }
    rows = []
    for name, term in terms.items():
        rows.append({"target": "CDK2", "prior_definition": name, **evaluate(raw, term)})
    return pd.DataFrame(rows)


def main() -> int:
    if not (RERUN / "FINGERPRINTS_DONE").exists():
        raise RuntimeError("Corrected 1PXN strict fingerprints are not complete")
    OUT.mkdir(parents=True, exist_ok=True)
    native, union, counts = primary_native_union()
    master = build_master(union)
    master.to_parquet(OUT / "cdk2_master.parquet", index=False)

    status = master.groupby(["target_receptor_id", "status"], dropna=False).size().rename("n_rows").reset_index()
    status.to_csv(OUT / "pose_fingerprint_status.csv", index=False)
    audit = pd.DataFrame([{
        "target": "CDK2",
        "primary_receptors_json": json.dumps(PRIMARY),
        "primary_native_receptors_json": json.dumps(PRIMARY),
        "expected_ligand_receptor_pairs": len(master),
        "n_ligands": master["lit_pcba_id"].nunique(),
        "strict_fingerprint_success": int(master["status"].eq("ok").sum()),
        "no_scored_pose": int(master["status"].eq("no_scored_pose").sum()),
        "strict_fingerprint_failure": int((~master["status"].isin(["ok", "no_scored_pose"])).sum()),
        "native_union_size": len(union),
    }])
    audit.to_csv(OUT / "cdk2_preparation_audit.csv", index=False)

    config = {
        "master": OUT / "cdk2_master.parquet",
        "actives": Path("/mnt/e/cdk2_dude/actives_final.ism"),
        "decoys": Path("/mnt/e/cdk2_dude/decoys_final.ism"),
        "target_size": len(union),
    }
    raw, per = statistics.load_target(config)
    summary, draws = statistics.run_nulls("CDK2", raw, config)
    intervals, effects = statistics.bootstrap_metric_tables("CDK2", per)
    summary.to_csv(OUT / "permutation_null_summary.csv", index=False)
    draws.to_parquet(OUT / "permutation_null_draws.parquet", index=False)
    intervals.to_csv(OUT / "bootstrap_metric_intervals.csv", index=False)
    effects.to_csv(OUT / "paired_metric_effects.csv", index=False)
    statistics.leave_one_receptor_out("CDK2", raw).to_csv(OUT / "leave_one_receptor_out.csv", index=False)
    statistics.formula_sensitivity("CDK2", raw).to_csv(OUT / "interaction_formula_sensitivity.csv", index=False)
    statistics.active_counts("CDK2", per).to_csv(OUT / "top1_active_counts.csv", index=False)
    prior_sensitivity(master, native, union, counts).to_csv(OUT / "native_prior_sensitivity.csv", index=False)
    bit_rows = [{
        "target": "CDK2", "interaction_bit": bit, "native_count": count,
        "native_frequency": count / len(native), "in_primary_union": True,
        "in_60pct_core": count / len(native) >= 0.60,
        "native_receptors_json": json.dumps(sorted(receptor for receptor, values in zip(native["receptor_id"], native["native_bits"]) if bit in values)),
    } for bit, count in sorted(counts.items())]
    pd.DataFrame(bit_rows).to_csv(OUT / "native_interaction_bits.csv", index=False)
    (OUT / "metadata.json").write_text(json.dumps({
        "source": "targeted_1pxn_altloc_rerun", "seed": statistics.SEED,
        "n_bootstrap": statistics.N_BOOT, "n_permutations": statistics.N_PERM,
        "primary_receptors": PRIMARY, "primary_native_receptors": PRIMARY,
        "native_union_size": len(union),
    }, indent=2) + "\n")
    (OUT / "REBUILD_DONE").write_text("done\n")
    print(summary.to_string(index=False), flush=True)
    print(effects.to_string(index=False), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
