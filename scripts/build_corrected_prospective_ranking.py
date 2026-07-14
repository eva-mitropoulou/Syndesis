"""Build the prospective ranking from strict, pose-coupled ProLIF fingerprints."""
from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from syndesis.enrichment.native_prior import native_union, parse_fingerprint  # noqa: E402

INPUTS = ROOT / "results/analysis_inputs"
FINGERPRINTS = INPUTS / "prospective_pose_fingerprints.parquet"
MASTER = INPUTS / "prospective_master_scores.parquet"
LIBRARY = INPUTS / "prospective_library.parquet"
OUTPUT = INPUTS / "prospective_ranked_corrected.parquet"


def main() -> int:
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    fingerprints = pd.read_parquet(FINGERPRINTS)
    failures = fingerprints[~fingerprints["status"].eq("ok")]
    if not failures.empty:
        raise RuntimeError(f"Prospective ranking blocked by {len(failures)} strict ProLIF failures")

    native = pd.read_parquet(INPUTS / "egfr_native_interaction_fingerprints.parquet")
    target, included_receptors = native_union(native, ["6duk_c_jbj_1103"])
    fingerprints["native_union_recall"] = fingerprints["fingerprint_sparse_json"].map(
        lambda value: len(parse_fingerprint(value) & target) / len(target)
    )
    fingerprints["zinc_id"] = fingerprints["ligand_id"].astype(str)

    master = pd.read_parquet(MASTER)
    master["zinc_id"] = master["lit_pcba_id"].astype(str)
    merged = master.merge(
        fingerprints[["zinc_id", "target_receptor_id", "native_union_recall", "num_interactions"]],
        on=["zinc_id", "target_receptor_id"], how="left", validate="one_to_one",
    )
    if merged["native_union_recall"].isna().any():
        raise RuntimeError(f"{int(merged['native_union_recall'].isna().sum())} prospective score rows lack fingerprints")
    merged["gnina_interaction_score"] = merged["cnnscore"] * (1.0 + merged["native_union_recall"])
    best = merged.sort_values("gnina_interaction_score", ascending=False).drop_duplicates("zinc_id")
    best = best.rename(columns={"native_union_recall": "key_interaction_recall_consensus"})

    library = pd.read_parquet(LIBRARY)
    library["zinc_id"] = library["zinc_id"].astype(str)
    ranked = best.merge(library, on="zinc_id", how="left", validate="one_to_one")
    ranked = ranked.sort_values("gnina_interaction_score", ascending=False).reset_index(drop=True)
    ranked["cnn_top_decile_threshold"] = ranked["cnnscore"].quantile(0.9)
    ranked["recall_median_threshold"] = ranked["key_interaction_recall_consensus"].median()
    ranked["passes_gate"] = (
        (ranked["cnnscore"] >= ranked["cnn_top_decile_threshold"])
        & (ranked["key_interaction_recall_consensus"] >= ranked["recall_median_threshold"])
    )
    ranked.to_parquet(OUTPUT, index=False)
    ranked.to_csv(OUTPUT.with_suffix(".csv"), index=False)
    print({
        "n_ranked": len(ranked),
        "n_gate_pass": int(ranked["passes_gate"].sum()),
        "cnn_threshold": float(ranked["cnn_top_decile_threshold"].iloc[0]),
        "recall_threshold": float(ranked["recall_median_threshold"].iloc[0]),
        "top_zinc_id": ranked.iloc[0]["zinc_id"],
        "native_prior_receptors": included_receptors,
        "native_prior_bit_count": len(target),
    })
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
