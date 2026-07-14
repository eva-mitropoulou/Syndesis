"""Measure EGFR active recovery by similarity to ATP-site native ligands."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem, DataStructs
from rdkit.Chem import rdFingerprintGenerator

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from syndesis.enrichment.metrics import enrichment_factor  # noqa: E402
from syndesis.enrichment.native_prior import native_union, parse_fingerprint, recall  # noqa: E402

INPUTS = ROOT / "results" / "analysis_inputs"
OUTPUT = ROOT / "results" / "robustness"
PRIMARY_EGFR_RECEPTORS = {
    "1m17_a_aq4_999", "1xkk_a_fmm_91", "4hjo_a_aq4_1001", "5cav_a_4zq_1101",
}


def rank_series(scores: pd.Series) -> pd.Series:
    return scores.rank(method="first", ascending=False).astype(int)


def main() -> int:
    OUTPUT.mkdir(parents=True, exist_ok=True)
    master = pd.read_parquet(INPUTS / "egfr_master.parquet")
    master = master[master["target_receptor_id"].isin(PRIMARY_EGFR_RECEPTORS)].copy()
    if set(master["target_receptor_id"].unique()) != PRIMARY_EGFR_RECEPTORS:
        raise RuntimeError("EGFR master scores do not contain exactly the four primary receptors")
    native = pd.read_parquet(INPUTS / "egfr_native_interaction_fingerprints.parquet")
    prior, _ = native_union(native, ["6duk_c_jbj_1103"])
    master["recall"] = master["fingerprint_sparse_json"].map(
        lambda value: recall(parse_fingerprint(value), prior)
    )
    grouped = master.groupby("lit_pcba_id", sort=True)
    ligands = pd.DataFrame({
        "label": grouped["label"].first(),
        "gnina_score": grouped["cnnscore"].max(),
        "coupled_score": master.assign(score=master["cnnscore"] * (1 + master["recall"])).groupby(
            "lit_pcba_id", sort=True
        )["score"].max(),
    }).reset_index()

    descriptors = pd.read_csv(INPUTS / "egfr_ligand_descriptors.csv", dtype={"lit_pcba_id": str})
    if "smiles" not in descriptors:
        raise RuntimeError("EGFR descriptor input must include canonical benchmark SMILES")
    ligands["lit_pcba_id"] = ligands["lit_pcba_id"].astype(str)
    ligands = ligands.merge(descriptors[["lit_pcba_id", "smiles"]], on="lit_pcba_id", validate="one_to_one")

    native_ligands = pd.read_csv(INPUTS / "egfr_native_ligands.csv")
    unique_native = native_ligands[native_ligands["primary_prior"]].drop_duplicates("inchi_key")
    generator = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    native_fps = [generator.GetFingerprint(Chem.MolFromSmiles(value)) for value in unique_native["ligand_smiles"]]

    active = ligands[ligands["label"].eq(1)].copy()
    active["max_native_ecfp4_tanimoto"] = [
        max(DataStructs.TanimotoSimilarity(generator.GetFingerprint(Chem.MolFromSmiles(smiles)), fp) for fp in native_fps)
        for smiles in active["smiles"]
    ]
    active["similarity_stratum"] = pd.cut(
        active["max_native_ecfp4_tanimoto"],
        bins=[-np.inf, 0.30, 0.50, 0.70, np.inf],
        labels=["<0.30", "0.30-<0.50", "0.50-<0.70", ">=0.70"],
        right=False,
    ).astype(str)
    ligands["gnina_rank"] = rank_series(ligands["gnina_score"])
    ligands["coupled_rank"] = rank_series(ligands["coupled_score"])
    active = active.merge(
        ligands[["lit_pcba_id", "gnina_rank", "coupled_rank"]], on="lit_pcba_id", validate="one_to_one"
    )
    top_n = max(1, int(round(len(ligands) * 0.01)))
    active["recovered_gnina_top1pct"] = active["gnina_rank"].le(top_n)
    active["recovered_coupled_top1pct"] = active["coupled_rank"].le(top_n)
    active.to_csv(OUTPUT / "egfr_active_native_similarity.csv", index=False)

    rows = []
    for stratum, group in active.groupby("similarity_stratum", observed=False):
        selected_ids = set(group["lit_pcba_id"])
        subset = ligands[ligands["label"].eq(0) | ligands["lit_pcba_id"].isin(selected_ids)]
        rows.append({
            "similarity_stratum": stratum,
            "n_actives": len(group),
            "mean_max_native_ecfp4_tanimoto": group["max_native_ecfp4_tanimoto"].mean(),
            "gnina_actives_in_global_top1pct": int(group["recovered_gnina_top1pct"].sum()),
            "coupled_actives_in_global_top1pct": int(group["recovered_coupled_top1pct"].sum()),
            "global_top1pct_recovery_difference": int(
                group["recovered_coupled_top1pct"].sum() - group["recovered_gnina_top1pct"].sum()
            ),
            "gnina_ef1_with_all_decoys": enrichment_factor(
                subset["label"].to_numpy(int), subset["gnina_score"].to_numpy(float), 0.01
            ),
            "coupled_ef1_with_all_decoys": enrichment_factor(
                subset["label"].to_numpy(int), subset["coupled_score"].to_numpy(float), 0.01
            ),
        })
    pd.DataFrame(rows).to_csv(OUTPUT / "egfr_active_similarity_strata.csv", index=False)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
