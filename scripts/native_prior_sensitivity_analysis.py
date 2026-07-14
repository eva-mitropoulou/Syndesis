"""Evaluate alternative native-derived interaction priors from archived pose bitsets."""
from __future__ import annotations

import argparse
import json
import sys
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from egfr_dockingforge.enrichment.metrics import bedroc, enrichment_factor, roc_auc  # noqa: E402


TARGETS = {
    "EGFR": {
        "master": ROOT / "results/analysis_inputs/egfr_master.parquet",
        "native": ROOT / "results/analysis_inputs/egfr_native_interaction_fingerprints.parquet",
        "exact_overlap_receptors": ["1xkk_a_fmm_91"],
        "descriptors": ROOT / "results/analysis_inputs/egfr_ligand_descriptors.csv",
    },
    "CDK2": {
        "master": ROOT / "results/analysis_inputs/cdk2_master.parquet",
        "native": ROOT / "results/analysis_inputs/cdk2_native_interaction_fingerprints.parquet",
        "exact_overlap_receptors": ["2a4l_a_rrc", "1aq1_a_stu"],
        "descriptors": ROOT / "results/analysis_inputs/cdk2_ligand_descriptors.csv",
    },
}


def bits(value: str) -> set[str]:
    parsed = json.loads(value) if isinstance(value, str) and value else []
    return {str(item) for item in parsed}


def recall(observed: set[str], target: set[str]) -> float:
    return len(observed & target) / len(target) if target else 0.0


def jaccard(observed: set[str], target: set[str]) -> float:
    union = observed | target
    return len(observed & target) / len(union) if union else 0.0


def weighted_recall(observed: set[str], weights: dict[str, float]) -> float:
    denominator = sum(weights.values())
    return sum(weight for bit, weight in weights.items() if bit in observed) / denominator if denominator else 0.0


def evaluate(raw: pd.DataFrame, term: pd.Series) -> dict[str, float]:
    scored = raw.assign(term=term, score=raw["cnnscore"] * (1.0 + term))
    grouped = scored.groupby("lit_pcba_id", sort=True)
    labels = grouped["label"].first().to_numpy(int)
    scores = grouped["score"].max().to_numpy(float)
    return {
        "roc_auc": roc_auc(labels, scores),
        "ef1": enrichment_factor(labels, scores, 0.01),
        "ef5": enrichment_factor(labels, scores, 0.05),
        "bedroc": bedroc(labels, scores, 80.5),
    }


def analyze_target(name: str, config: dict, out: Path) -> None:
    native = pd.read_parquet(config["native"]).copy()
    native["native_bits"] = native["fingerprint_sparse_json"].map(bits)
    native_sets = dict(zip(native["receptor_id"], native["native_bits"]))
    counts = Counter(bit for value in native["native_bits"] for bit in value)
    union = set(counts)
    core = {bit for bit, count in counts.items() if count / len(native) >= 0.60}
    weights = {bit: count / len(native) for bit, count in counts.items()}

    bit_rows = []
    for bit in sorted(union):
        bit_rows.append({
            "target": name,
            "interaction_bit": bit,
            "native_count": counts[bit],
            "native_frequency": weights[bit],
            "in_global_union": True,
            "in_60pct_core": bit in core,
            "native_receptors_json": json.dumps(sorted(r for r, value in native_sets.items() if bit in value)),
        })
    pd.DataFrame(bit_rows).to_csv(out / f"{name.lower()}_native_interaction_bits.csv", index=False)

    merged = pd.read_parquet(config["master"]).sort_values("cnnscore", ascending=False).drop_duplicates(
        ["lit_pcba_id", "target_receptor_id"]
    )
    merged["lit_pcba_id"] = merged["lit_pcba_id"].astype(str)
    allowed = merged["status"].isin(["ok", "no_scored_pose"])
    if not allowed.all():
        raise RuntimeError(f"{name}: {int((~allowed).sum())} strict fingerprints failed for scored poses")
    ok = merged["status"].eq("ok")
    merged["pose_bits"] = None
    merged.loc[ok, "pose_bits"] = merged.loc[ok, "fingerprint_sparse_json"].map(bits)
    failure_report = (
        merged.assign(error_type=merged["error"].fillna("").str.split(":", n=1).str[0])
        .groupby(["target_receptor_id", "status", "error_type"], dropna=False)
        .size().rename("n_rows").reset_index()
    )
    failure_report.insert(0, "target", name)
    failure_report.to_csv(out / f"{name.lower()}_pose_fingerprint_status.csv", index=False)

    terms: dict[str, pd.Series] = {
        "global_union_recall": merged["pose_bits"].map(lambda value: recall(value, union) if isinstance(value, set) else np.nan),
        "60pct_core_recall": merged["pose_bits"].map(lambda value: recall(value, core) if isinstance(value, set) else np.nan),
        "frequency_weighted_recall": merged["pose_bits"].map(lambda value: weighted_recall(value, weights) if isinstance(value, set) else np.nan),
        "global_union_jaccard": merged["pose_bits"].map(lambda value: jaccard(value, union) if isinstance(value, set) else np.nan),
        "receptor_specific_recall": pd.Series(
            [recall(value, native_sets[receptor]) if isinstance(value, set) else np.nan for value, receptor in zip(merged["pose_bits"], merged["target_receptor_id"])],
            index=merged.index,
        ),
    }

    descriptors = pd.read_csv(config["descriptors"], dtype={"lit_pcba_id": str}).set_index("lit_pcba_id")
    merged["heavy_atom_count"] = merged["lit_pcba_id"].map(descriptors["heavy_atom_count"])
    merged["molecular_weight"] = merged["lit_pcba_id"].map(descriptors["molecular_weight"])
    if merged[["heavy_atom_count", "molecular_weight"]].isna().any().any():
        raise RuntimeError(f"{name}: ligand descriptor table is incomplete")
    merged["global_union_recall"] = terms["global_union_recall"]
    correlation_rows = []
    for scope, frame in {
        "all_receptor_poses": merged,
        "best_coupled_pose_per_ligand": merged.assign(
            coupled=merged["cnnscore"] * (1 + merged["global_union_recall"])
        ).sort_values("coupled", ascending=False).drop_duplicates("lit_pcba_id"),
    }.items():
        for variable in ["heavy_atom_count", "molecular_weight", "num_interactions"]:
            correlation_rows.append({
                "target": name,
                "scope": scope,
                "interaction_term": "global_union_recall",
                "variable": variable,
                "spearman_rho": frame["global_union_recall"].corr(frame[variable], method="spearman"),
                "n_rows": len(frame),
            })
    correlation_path = out / "interaction_size_correlations.csv"
    correlations = pd.DataFrame(correlation_rows)
    if correlation_path.exists():
        correlations = pd.concat([pd.read_csv(correlation_path), correlations], ignore_index=True)
        correlations = correlations.drop_duplicates(["target", "scope", "interaction_term", "variable"], keep="last")
    correlations.to_csv(correlation_path, index=False)

    rows = []
    for definition, term in terms.items():
        rows.append({"target": name, "prior_definition": definition, "excluded_native_receptor": "", **evaluate(merged, term)})

    for excluded, excluded_bits in native_sets.items():
        loo_union = set().union(*(value for receptor, value in native_sets.items() if receptor != excluded))
        term = merged["pose_bits"].map(lambda value, target=loo_union: recall(value, target) if isinstance(value, set) else np.nan)
        rows.append({
            "target": name,
            "prior_definition": "leave_one_native_out_union_recall",
            "excluded_native_receptor": excluded,
            **evaluate(merged, term),
        })

    exact = config["exact_overlap_receptors"]
    if exact:
        excluded_union = set().union(*(value for receptor, value in native_sets.items() if receptor not in exact))
        term = merged["pose_bits"].map(lambda value: recall(value, excluded_union) if isinstance(value, set) else np.nan)
        rows.append({
            "target": name,
            "prior_definition": "all_exact_overlap_natives_excluded_union_recall",
            "excluded_native_receptor": ";".join(exact),
            **evaluate(merged, term),
        })

    pd.DataFrame(rows).to_csv(out / f"{name.lower()}_native_prior_sensitivity.csv", index=False)
    summary = {
        "target": name,
        "n_native_complexes": len(native),
        "global_union_size": len(union),
        "core_60pct_size": len(core),
        "n_pose_rows": len(merged),
        "n_fingerprint_success": int(ok.sum()),
        "n_no_scored_pose": int(merged["status"].eq("no_scored_pose").sum()),
        "n_fingerprint_failure": int((~allowed).sum()),
        "failure_policy": "Fail closed for every scored pose; rows without a docked/scored pose remain explicitly non-applicable and are not imputed.",
    }
    (out / f"{name.lower()}_native_prior_metadata.json").write_text(json.dumps(summary, indent=2) + "\n")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-root", type=Path, default=ROOT / "results/robustness")
    parser.add_argument("--target", choices=sorted(TARGETS), action="append")
    args = parser.parse_args()
    args.output_root.mkdir(parents=True, exist_ok=True)
    selected = args.target or list(TARGETS)
    for name in selected:
        config = TARGETS[name]
        analyze_target(name, config, args.output_root)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
