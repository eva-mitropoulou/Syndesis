from __future__ import annotations

import json
import math
from collections import Counter
from typing import Any

import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Crippen, Descriptors, Lipinski, rdMolDescriptors

from egfr_dockingforge.common.io import write_json, write_table
from egfr_dockingforge.stage6 import schemas


def _parse_bits(value: Any) -> list[str]:
    if value is None or (isinstance(value, float) and math.isnan(value)):
        return []
    if isinstance(value, list):
        return [str(x) for x in value]
    text = str(value)
    if not text:
        return []
    try:
        payload = json.loads(text)
        if isinstance(payload, list):
            return [str(x) for x in payload]
    except json.JSONDecodeError:
        pass
    return [x for x in text.split("|") if x]


def _ligand_descriptors(ligands: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for row in ligands.to_dict("records"):
        mol = None
        for key in ("prepared_ligand_file", "native_ligand_file", "immutable_reference_pose_file"):
            path = row.get(key)
            if isinstance(path, str) and path:
                mol = Chem.MolFromPDBFile(path, removeHs=False, sanitize=True)
                if mol is not None:
                    break
        if mol is None:
            raise RuntimeError(f"RDKit could not read a ligand structure for {row.get('ligand_id')}")
        rows.append(
            {
                "ligand_id": row["ligand_id"],
                "ligand_heavy_atom_count": int(mol.GetNumHeavyAtoms()),
                "ligand_rotatable_bonds": int(Lipinski.NumRotatableBonds(mol)),
                "ligand_formal_charge": int(sum(atom.GetFormalCharge() for atom in mol.GetAtoms())),
                "ligand_aromatic_ring_count": int(rdMolDescriptors.CalcNumAromaticRings(mol)),
                "ligand_molecular_weight": float(Descriptors.MolWt(mol)),
                "ligand_tpsa": float(rdMolDescriptors.CalcTPSA(mol)),
                "ligand_logp": float(Crippen.MolLogP(mol)),
                "ligand_hbd": int(Lipinski.NumHDonors(mol)),
                "ligand_hba": int(Lipinski.NumHAcceptors(mol)),
                "ligand_source": str(row.get("ligand_prep_tool", "stage3")),
            }
        )
    return pd.DataFrame(rows)


def _fingerprint_bit_columns(fps: pd.DataFrame, min_frequency: int) -> pd.DataFrame:
    bit_lists = fps["fingerprint_sparse_json"].map(_parse_bits)
    counts = Counter(bit for bits in bit_lists for bit in bits)
    selected = sorted(bit for bit, count in counts.items() if count >= min_frequency)
    frame = pd.DataFrame({"pose_id": fps["pose_id"]})
    for bit in selected:
        name = "ifp_bit_" + "".join(ch if ch.isalnum() else "_" for ch in bit).strip("_").lower()
        frame[name] = bit_lists.map(lambda bits, b=bit: int(b in bits))
    return frame


def build_pose_features(inputs: dict[str, pd.DataFrame], config: dict[str, Any], paths: dict[str, Any]) -> pd.DataFrame:
    base = inputs["stage5_interaction_features"].copy()
    score = inputs["pose_score_table"][
        [
            "pose_id",
            "vina_rescore",
            "vinardo_rescore",
            "gnina_empirical_affinity",
            "cnnscore",
            "cnnaffinity",
            "cnn_vs",
        ]
    ].copy()
    sanity = inputs["pose_sanity"].copy()
    recovery = inputs["interaction_recovery"][
        [
            "pose_id",
            "key_interaction_recall_consensus",
            "key_interaction_precision_consensus",
            "key_interaction_f1_consensus",
            "ifp_tanimoto_to_consensus",
            "hinge_interaction_recovered_flag",
            "catalytic_lys_glu_region_consistent_flag",
            "gatekeeper_region_consistent_flag",
            "dfg_region_consistent_flag",
            "interaction_recovery_label",
        ]
    ].copy()
    clusters = inputs["binding_mode_clusters"]
    cluster_cols = clusters.loc[clusters["entity_type"].eq("docked_pose"), ["pose_id", "distance_to_medoid", "cluster_label"]].copy()
    receptors = inputs["receptor_ensemble"][
        [
            "receptor_id",
            "dfg_state",
            "chelix_state",
            "saltbridge_state",
            "hrd_state",
            "activation_loop_state",
            "mutation_flag",
            "active_site_completeness_score",
            "cluster_id",
            "suggested_docking_box_size",
        ]
    ].rename(columns={"receptor_id": "target_receptor_id", "cluster_id": "receptor_cluster_id"})
    ligands = _ligand_descriptors(inputs["ligand_prep"])
    fp_bits = _fingerprint_bit_columns(inputs["docked_pose_fingerprints"], int(config["features"]["min_fingerprint_bit_frequency"]))

    df = base.merge(score, on="pose_id", how="left", suffixes=("", "_stage4"))
    df = df.merge(sanity, on="pose_id", how="left", suffixes=("", "_sanity"))
    df = df.merge(recovery, on="pose_id", how="left", suffixes=("", "_recovery"))
    df = df.merge(cluster_cols, on="pose_id", how="left")
    df = df.merge(receptors, on="target_receptor_id", how="left")
    df = df.merge(ligands, on="ligand_id", how="left")
    df = df.merge(fp_bits, on="pose_id", how="left")

    df["native_receptor_id"] = df["ligand_id"].str.replace(r"_[^_]+$", "", regex=True)
    df["pose_rank"] = df["pose_rank"].astype(int)
    df["score_disagreement_cnn_minus_docking"] = pd.to_numeric(df["gnina_cnnscore"], errors="coerce") - pd.to_numeric(df["original_docking_score"], errors="coerce")
    df["score_disagreement_affinity_minus_docking"] = pd.to_numeric(df["gnina_cnnaffinity"], errors="coerce") - pd.to_numeric(df["original_docking_score"], errors="coerce")
    df["rank_fraction_within_task"] = df.groupby("docking_task_id")["pose_rank"].rank(method="first") / df.groupby("docking_task_id")["pose_rank"].transform("count")
    df["sanity_status_encoded"] = df["sanity_status"].map({"pass": 1, "warning": 0.5, "fail": 0}).fillna(0)
    df["dominant_binding_mode_cluster_compatibility"] = df["cluster_label"].eq("native_enriched").astype(int)
    df["binding_mode_cluster_distance"] = pd.to_numeric(df["distance_to_medoid"], errors="coerce")
    df["state_match_flag"] = [
        int("native" in str(task_type) and str(target).split("_")[0] == str(ligand).split("_")[0])
        for task_type, target, ligand in zip(df["task_type"], df["target_receptor_id"], df["ligand_id"], strict=False)
    ]
    size = df["suggested_docking_box_size"].fillna("[0,0,0]").astype(str).str.extract(r"\[?([0-9.]+),\s*([0-9.]+),\s*([0-9.]+)\]?")
    df["docking_box_volume"] = size.astype(float).prod(axis=1)
    df["mutation_flag"] = df["mutation_flag"].fillna(False).astype(bool)

    if "gnina_cnnscore" in df and "cnnscore" in df:
        df["cnnscore"] = df["cnnscore"].fillna(df["gnina_cnnscore"])
    if "gnina_cnnaffinity" in df and "cnnaffinity" in df:
        df["cnnaffinity"] = df["cnnaffinity"].fillna(df["gnina_cnnaffinity"])

    for col in df.columns:
        if df[col].dtype == bool:
            df[col] = df[col].astype(int)
    df = df.loc[:, ~df.columns.duplicated()].copy()
    write_table(paths["processed"] / "pose_model_features.parquet", df)
    write_table(paths["processed"] / "pose_model_features.csv", df)
    write_json(paths["processed"] / "pose_model_feature_summary.json", {"rows": int(len(df)), "columns": int(len(df.columns))})
    return df


def training_feature_columns(features: pd.DataFrame, audit: pd.DataFrame) -> list[str]:
    metadata = set(schemas.IDENTIFIER_COLUMNS) | {"fingerprint_sparse_json", "fingerprint_bitstring", "warnings_json", "parse_warnings_json", "sanity_warnings_json"}
    allowed = set(audit.loc[audit["allowed_for_training"], "feature_name"])
    columns: list[str] = []
    for col in features.columns:
        if col in metadata or col not in allowed:
            continue
        if features[col].dtype == object:
            if col.endswith("_json"):
                continue
        columns.append(col)
    return columns
