"""Generate submission-facing robustness tables from archived enrichment results."""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from syndesis.enrichment.metrics import bedroc, enrichment_factor, roc_auc  # noqa: E402
from syndesis.enrichment.native_prior import (  # noqa: E402
    jaccard,
    native_union,
    parse_fingerprint,
    recall,
)

OUT = ROOT / "results" / "robustness"
OUT.mkdir(parents=True, exist_ok=True)
SEED = 807
N_PERM = 1000
N_BOOT = 2000

TARGETS = {
    "EGFR": {
        "master": ROOT / "results/analysis_inputs/egfr_master.parquet",
        "native": ROOT / "results/analysis_inputs/egfr_native_interaction_fingerprints.parquet",
        "descriptors": ROOT / "results/analysis_inputs/egfr_ligand_descriptors.csv",
        "excluded_native_receptors": ["6duk_c_jbj_1103"],
        "docking_receptors": [
            "1m17_a_aq4_999", "1xkk_a_fmm_91", "4hjo_a_aq4_1001", "5cav_a_4zq_1101",
        ],
        "analysis_role": "primary_four_receptor",
    },
    "EGFR_5RECEPTOR_SENSITIVITY": {
        "master": ROOT / "results/analysis_inputs/egfr_master.parquet",
        "native": ROOT / "results/analysis_inputs/egfr_native_interaction_fingerprints.parquet",
        "descriptors": ROOT / "results/analysis_inputs/egfr_ligand_descriptors.csv",
        "excluded_native_receptors": ["6duk_c_jbj_1103"],
        "docking_receptors": [
            "1m17_a_aq4_999", "1xkk_a_fmm_91", "4hjo_a_aq4_1001", "5cav_a_4zq_1101",
            "6duk_c_jbj_1103",
        ],
        "analysis_role": "five_receptor_sensitivity",
    },
    "CDK2_5RECEPTOR_SENSITIVITY": {
        "master": ROOT / "results/analysis_inputs/cdk2_master.parquet",
        "native": ROOT / "results/analysis_inputs/cdk2_native_interaction_fingerprints.parquet",
        "descriptors": ROOT / "results/analysis_inputs/cdk2_ligand_descriptors.csv",
        "excluded_native_receptors": [],
        "docking_receptors": [
            "1qmz_a_atp", "1fin_a_atp", "2a4l_a_rrc", "1aq1_a_stu", "1pxn_a_ck6",
        ],
        "analysis_role": "five_receptor_sensitivity_unphosphorylated_1qmz",
    },
    "CDK2": {
        "master": ROOT / "results/analysis_inputs/cdk2_master.parquet",
        "native": ROOT / "results/analysis_inputs/cdk2_native_interaction_fingerprints.parquet",
        "descriptors": ROOT / "results/analysis_inputs/cdk2_ligand_descriptors.csv",
        "excluded_native_receptors": [],
        "docking_receptors": [
            "1fin_a_atp", "2a4l_a_rrc", "1aq1_a_stu", "1pxn_a_ck6",
        ],
        "analysis_role": "primary_four_receptor_excluding_unphosphorylated_1qmz",
    },
}


def ef(labels: np.ndarray, scores: np.ndarray, fraction: float = 0.01) -> float:
    return enrichment_factor(labels, scores, fraction)


def load_target(config: dict) -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    raw = pd.read_parquet(config["master"])
    expected_receptors = set(config["docking_receptors"])
    raw = raw[raw["target_receptor_id"].isin(expected_receptors)].copy()
    observed_receptors = set(raw["target_receptor_id"].unique())
    if observed_receptors != expected_receptors:
        raise RuntimeError(
            f"Receptor filter mismatch: expected {sorted(expected_receptors)}, "
            f"observed {sorted(observed_receptors)}"
        )
    raw = raw.sort_values("cnnscore", ascending=False).drop_duplicates(
        ["lit_pcba_id", "target_receptor_id"]
    )
    allowed = raw["status"].isin(["ok", "no_scored_pose"])
    if not allowed.all():
        raise RuntimeError(f"{int((~allowed).sum())} scored poses have failed fingerprints")
    native = pd.read_parquet(config["native"])
    target, included = native_union(native, config["excluded_native_receptors"])
    raw["pose_bits"] = None
    scored = raw["status"].eq("ok")
    raw.loc[scored, "pose_bits"] = raw.loc[scored, "fingerprint_sparse_json"].map(parse_fingerprint)
    raw["recall"] = raw["pose_bits"].map(
        lambda value: recall(value, target) if isinstance(value, set) else np.nan
    )
    raw["tanimoto"] = raw["pose_bits"].map(
        lambda value: jaccard(value, target) if isinstance(value, set) else np.nan
    )

    if "num_interactions" not in raw:
        raise RuntimeError("Corrected enrichment master must contain num_interactions")

    intersection = raw["pose_bits"].map(
        lambda value: len(value & target) if isinstance(value, set) else np.nan
    )
    raw["precision"] = np.where(raw["num_interactions"] > 0, intersection / raw["num_interactions"], 0.0)
    denom = raw["precision"] + raw["recall"]
    raw["f1"] = np.where(denom > 0, 2 * raw["precision"] * raw["recall"] / denom, 0.0)

    grouped = raw.groupby("lit_pcba_id", sort=True)
    per = pd.DataFrame({"label": grouped["label"].first()})
    per["gnina"] = grouped["cnnscore"].max()
    per["coupled"] = raw.assign(score=raw["cnnscore"] * (1 + raw["recall"])).groupby(
        "lit_pcba_id"
    )["score"].max()
    prior = {
        "included_native_receptors": included,
        "excluded_native_receptors": config["excluded_native_receptors"],
        "interaction_bit_count": len(target),
        "docking_receptors": sorted(expected_receptors),
        "analysis_role": config["analysis_role"],
    }
    return raw, per.reset_index(), prior


def aligned_vectors(raw: pd.DataFrame):
    cnn_frame = raw.pivot(index="lit_pcba_id", columns="target_receptor_id", values="cnnscore").sort_index()
    recall_frame = raw.pivot(index="lit_pcba_id", columns="target_receptor_id", values="recall").reindex_like(cnn_frame)
    ligand_ids = cnn_frame.index.astype(str).tolist()
    cnn = cnn_frame.to_numpy(float)
    recall = recall_frame.to_numpy(float)
    if cnn.shape[1] != raw["target_receptor_id"].nunique() or recall.shape != cnn.shape:
        raise RuntimeError("Every ligand must have one score and fingerprint for every receptor")
    labels = raw.groupby("lit_pcba_id", sort=True)["label"].first().reindex(cnn_frame.index).to_numpy(int)
    return ligand_ids, cnn, recall, labels


def run_nulls(name: str, raw: pd.DataFrame, config: dict) -> tuple[pd.DataFrame, pd.DataFrame]:
    ligand_ids, cnn, recall, labels = aligned_vectors(raw)
    order = np.arange(len(ligand_ids))
    observed_scores = np.nanmax(cnn * (1 + recall), axis=1)
    observed = ef(labels, observed_scores)

    descriptors = pd.read_csv(config["descriptors"], dtype={"lit_pcba_id": str}).set_index("lit_pcba_id")
    heavy_atoms = descriptors["heavy_atom_count"].reindex(ligand_ids).to_numpy(float)
    if np.isnan(heavy_atoms).any():
        raise RuntimeError("Ligand descriptor table does not cover every benchmark molecule")
    bins = np.full(len(ligand_ids), -1, dtype=int)
    valid = ~np.isnan(heavy_atoms)
    bins[valid] = pd.qcut(heavy_atoms[valid], 10, labels=False, duplicates="drop")

    complete = np.isfinite(cnn).all(axis=1) & np.isfinite(recall).all(axis=1)
    definitions = {
        "all_ligand": [np.where(complete)[0]],
        "heavy_atom_decile": [np.where(complete & (bins == b))[0] for b in np.unique(bins)],
        "class_conditional": [np.where(complete & (labels == label))[0] for label in np.unique(labels)],
    }
    summaries = []
    draws = []
    for null_name, strata in definitions.items():
        values = []
        for permutation in range(N_PERM):
            rng = np.random.default_rng(SEED + permutation)
            assignment = order.copy()
            for indices in strata:
                assignment[indices] = rng.permutation(indices)
            scores = np.nanmax(cnn * (1 + recall[assignment]), axis=1)
            value = ef(labels, scores)
            values.append(value)
            draws.append({"target": name, "null": null_name, "permutation": permutation, "ef1": value})
        array = np.asarray(values)
        summaries.append({
            "target": name,
            "null": null_name,
            "n_permutations": N_PERM,
            "observed_ef1": observed,
            "null_mean": array.mean(),
            "null_sd": array.std(),
            "null_p2_5": np.percentile(array, 2.5),
            "null_p97_5": np.percentile(array, 97.5),
            "observed_minus_null_mean": observed - array.mean(),
            "empirical_p": (np.sum(array >= observed) + 1) / (N_PERM + 1),
        })
    return pd.DataFrame(summaries), pd.DataFrame(draws)


def paired_bootstrap_delta(labels: np.ndarray, baseline: np.ndarray, combined: np.ndarray) -> tuple[float, float, float]:
    positive = np.flatnonzero(labels == 1)
    negative = np.flatnonzero(labels == 0)
    values = []
    for iteration in range(N_BOOT):
        rng = np.random.default_rng(SEED + iteration)
        indices = np.concatenate([
            rng.choice(positive, len(positive), replace=True),
            rng.choice(negative, len(negative), replace=True),
        ])
        values.append(ef(labels[indices], combined[indices]) - ef(labels[indices], baseline[indices]))
    values = np.asarray(values)
    return np.percentile(values, 2.5), np.percentile(values, 97.5), np.mean(values > 0)


def metric_values(labels: np.ndarray, scores: np.ndarray) -> dict[str, float]:
    return {
        "roc_auc": roc_auc(labels, scores),
        "ef1": ef(labels, scores),
        "ef5": ef(labels, scores, 0.05),
        "bedroc_80_5": bedroc(labels, scores, 80.5),
    }


def bootstrap_metric_tables(
    name: str, per: pd.DataFrame
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Return class-stratified intervals for each arm and paired arm differences."""
    labels = per["label"].to_numpy(int)
    scores = {arm: per[arm].to_numpy(float) for arm in ["gnina", "coupled"]}
    positive = np.flatnonzero(labels == 1)
    negative = np.flatnonzero(labels == 0)
    arm_draws = {arm: {metric: [] for metric in metric_values(labels, values)} for arm, values in scores.items()}
    delta_draws = {metric: [] for metric in metric_values(labels, scores["gnina"])}
    for iteration in range(N_BOOT):
        rng = np.random.default_rng(SEED + iteration)
        indices = np.concatenate([
            rng.choice(positive, len(positive), replace=True),
            rng.choice(negative, len(negative), replace=True),
        ])
        sampled = {
            arm: metric_values(labels[indices], values[indices])
            for arm, values in scores.items()
        }
        for arm in scores:
            for metric, value in sampled[arm].items():
                arm_draws[arm][metric].append(value)
        for metric in delta_draws:
            delta_draws[metric].append(sampled["coupled"][metric] - sampled["gnina"][metric])

    intervals = []
    effects = []
    points = {arm: metric_values(labels, values) for arm, values in scores.items()}
    for arm in scores:
        for metric, estimate in points[arm].items():
            values = np.asarray(arm_draws[arm][metric])
            intervals.append({
                "target": name,
                "arm": arm,
                "metric": metric,
                "estimate": estimate,
                "ci_lo": np.percentile(values, 2.5),
                "ci_hi": np.percentile(values, 97.5),
                "n_bootstrap": N_BOOT,
            })
    for metric, values_list in delta_draws.items():
        values = np.asarray(values_list)
        effects.append({
            "target": name,
            "contrast": "coupled_minus_gnina",
            "metric": metric,
            "estimate": points["coupled"][metric] - points["gnina"][metric],
            "ci_lo": np.percentile(values, 2.5),
            "ci_hi": np.percentile(values, 97.5),
            "bootstrap_fraction_positive": np.mean(values > 0),
            "n_bootstrap": N_BOOT,
        })
    return pd.DataFrame(intervals), pd.DataFrame(effects)


def leave_one_receptor_out(name: str, raw: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for excluded in sorted(raw["target_receptor_id"].unique()):
        subset = raw[raw["target_receptor_id"] != excluded].copy()
        grouped = subset.groupby("lit_pcba_id", sort=True)
        labels = grouped["label"].first().to_numpy(int)
        baseline = grouped["cnnscore"].max().to_numpy(float)
        combined = subset.assign(score=subset["cnnscore"] * (1 + subset["recall"])).groupby(
            "lit_pcba_id", sort=True
        )["score"].max().to_numpy(float)
        lo, hi, frac = paired_bootstrap_delta(labels, baseline, combined)
        rows.append({
            "target": name,
            "excluded_receptor": excluded,
            "n_ligands": len(labels),
            "gnina_ef1": ef(labels, baseline),
            "coupled_ef1": ef(labels, combined),
            "delta_ef1": ef(labels, combined) - ef(labels, baseline),
            "delta_ci_lo": lo,
            "delta_ci_hi": hi,
            "bootstrap_fraction_positive": frac,
        })
    return pd.DataFrame(rows)


def formula_sensitivity(name: str, raw: pd.DataFrame) -> pd.DataFrame:
    grouped = raw.groupby("lit_pcba_id", sort=True)
    labels = grouped["label"].first().to_numpy(int)
    rows = []
    for term in ["recall", "tanimoto", "precision", "f1"]:
        for weight in [0.0, 0.25, 0.5, 1.0, 1.5, 2.0, 3.0]:
            cnn_rank = raw["cnnscore"].rank(method="average", pct=True)
            term_rank = raw[term].rank(method="average", pct=True)
            combinations = {
                "multiplicative": raw["cnnscore"] * (1 + weight * raw[term]),
                "additive": raw["cnnscore"] + weight * raw[term],
                "rank_fusion": cnn_rank + weight * term_rank,
            }
            for combination, pose_scores in combinations.items():
                scores = raw.assign(score=pose_scores).groupby(
                    "lit_pcba_id", sort=True
                )["score"].max().to_numpy(float)
                rows.append({
                    "target": name,
                    "combination": combination,
                    "term": term,
                    "lambda": weight,
                    "roc_auc": roc_auc(labels, scores),
                    "ef1": ef(labels, scores),
                    "ef5": ef(labels, scores, 0.05),
                    "bedroc": bedroc(labels, scores, 80.5),
                })
    return pd.DataFrame(rows)


def active_counts(name: str, per: pd.DataFrame) -> pd.DataFrame:
    labels = per["label"].to_numpy(int)
    n_top = max(1, int(round(len(per) * 0.01)))
    rows = []
    for arm in ["gnina", "coupled"]:
        order = np.argsort(-per[arm].to_numpy(float), kind="mergesort")
        rows.append({
            "target": name,
            "arm": arm,
            "n_ligands": len(per),
            "top_fraction": 0.01,
            "n_top": n_top,
            "n_actives_top": int(labels[order][:n_top].sum()),
            "n_actives_total": int(labels.sum()),
            "ef1": ef(labels, per[arm].to_numpy(float)),
        })
    return pd.DataFrame(rows)


def prospective_audit() -> pd.DataFrame:
    ranked = pd.read_parquet(ROOT / "results/analysis_inputs/prospective_ranked_corrected.parquet")
    cnn_threshold = ranked["cnn_top_decile_threshold"].iloc[0]
    recall_threshold = ranked["recall_median_threshold"].iloc[0]
    accepted = ranked[ranked["passes_gate"]].sort_values("gnina_interaction_score", ascending=False)
    top = ranked.sort_values("gnina_interaction_score", ascending=False).iloc[0]
    summary = {
        "n_scored": len(ranked),
        "cnn_top_decile_threshold": cnn_threshold,
        "global_median_recall_threshold": recall_threshold,
        "n_gate_pass": int(ranked["passes_gate"].sum()),
        "ranking_after_gate": "gnina_interaction_score_descending",
        "top_coupled_zinc_id": str(top["zinc_id"]),
        "top_coupled_cnnscore": top["cnnscore"],
        "top_coupled_recall": top["key_interaction_recall_consensus"],
        "top_coupled_score": top["gnina_interaction_score"],
        "top_coupled_passes_gate": bool(top["passes_gate"]),
        "top_accepted_zinc_id": str(accepted.iloc[0]["zinc_id"]),
    }
    accepted.to_csv(OUT / "prospective_gate_pass_candidates.csv", index=False)
    return pd.DataFrame([summary])


def analog_lineage() -> pd.DataFrame:
    candidates = pd.read_csv(ROOT / "data/processed/stage9/analog_candidates.csv")
    screening = pd.read_csv(ROOT / "data/processed/stage9/analog_screening_results.csv")
    acceptance = pd.read_csv(ROOT / "data/processed/stage9/analog_acceptance.csv")
    deterministic = candidates[candidates["generated_by"] == "rdkit_deterministic_transform"].copy()
    lineage = deterministic.merge(
        screening[["analog_id", "best_pose_id", "best_receptor_id"]], on="analog_id", how="left"
    ).merge(
        acceptance[["analog_id", "accepted_flag", "acceptance_tier", "score_hacking_flag", "rejection_reason", "acceptance_reason"]],
        on="analog_id",
        how="left",
    )
    columns = [
        "analog_id", "parent_molecule_id", "parent_smiles", "analog_smiles",
        "transformation_class", "generated_by", "source", "accepted_flag",
        "acceptance_tier", "score_hacking_flag", "rejection_reason", "acceptance_reason",
        "best_pose_id", "best_receptor_id",
    ]
    return lineage[columns]


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--target", choices=sorted(TARGETS), action="append")
    parser.add_argument(
        "--combine-only",
        action="store_true",
        help="Combine existing per-target artifacts without rerunning their resampling loops.",
    )
    args = parser.parse_args()
    if args.combine_only and args.target:
        parser.error("--combine-only cannot be combined with --target")
    selected = [] if args.combine_only else (args.target or list(TARGETS))
    output_names = {
        "permutation_null_summary": "csv",
        "permutation_null_draws": "parquet",
        "leave_one_receptor_out": "csv",
        "interaction_formula_sensitivity": "csv",
        "top1_active_counts": "csv",
        "bootstrap_metric_intervals": "csv",
        "paired_metric_effects": "csv",
    }
    for name in selected:
        config = TARGETS[name]
        raw, per, prior = load_target(config)
        summary, draws = run_nulls(name, raw, config)
        frames = {
            "permutation_null_summary": summary,
            "permutation_null_draws": draws,
            "leave_one_receptor_out": leave_one_receptor_out(name, raw),
            "interaction_formula_sensitivity": formula_sensitivity(name, raw),
            "top1_active_counts": active_counts(name, per),
        }
        intervals, effects = bootstrap_metric_tables(name, per)
        frames["bootstrap_metric_intervals"] = intervals
        frames["paired_metric_effects"] = effects
        prefix = f"interim_{name.lower()}_"
        for stem, frame in frames.items():
            path = OUT / f"{prefix}{stem}.{output_names[stem]}"
            if output_names[stem] == "parquet":
                frame.to_parquet(path, index=False)
            else:
                frame.to_csv(path, index=False)
        print(summary.to_string(index=False), flush=True)
        config["primary_prior_metadata"] = prior

    for stem, extension in output_names.items():
        paths = [OUT / f"interim_{name.lower()}_{stem}.{extension}" for name in TARGETS]
        if not all(path.exists() for path in paths):
            continue
        frames = [pd.read_parquet(path) if extension == "parquet" else pd.read_csv(path) for path in paths]
        combined = pd.concat(frames, ignore_index=True)
        if extension == "parquet":
            combined.to_parquet(OUT / f"{stem}.{extension}", index=False)
            if stem == "permutation_null_draws":
                combined.to_parquet(
                    ROOT / "results/analysis_inputs/permutation_null_draws.parquet",
                    index=False,
                )
        else:
            combined.to_csv(OUT / f"{stem}.{extension}", index=False)
    prospective_audit().to_csv(OUT / "prospective_gate_audit.csv", index=False)

    metadata = {
        "seed": SEED,
        "n_permutations": N_PERM,
        "n_bootstrap": N_BOOT,
        "ef_top_set_rounding": "round(N * fraction), minimum 1",
        "tie_handling": "stable mergesort descending",
        "primary_native_priors": {
            name: config.get("primary_prior_metadata") for name, config in TARGETS.items()
        },
        "analysis_roles": {name: config["analysis_role"] for name, config in TARGETS.items()},
        "ensemble_definitions": {
            name: {
                "docking_receptors": config["docking_receptors"],
                "excluded_native_receptors": config["excluded_native_receptors"],
            }
            for name, config in TARGETS.items()
        },
        "paper_analysis_config": "configs/paper_analysis.yaml",
    }
    (OUT / "analysis_metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")
    print(prospective_audit().to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
