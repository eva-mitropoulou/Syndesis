from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pandas as pd

from syndesis.common.io import write_table
from syndesis.stage5.prolif_engine import compute_interactions, fingerprint_from_interactions
from syndesis.stage5.residue_mapping import normalize_interaction_key
from syndesis.stage5.schemas import (
    KEY_INTERACTION_COLUMNS,
    NATIVE_FINGERPRINT_COLUMNS,
    NATIVE_INTERACTIONS_LONG_COLUMNS,
)


def _state_for_complex(ensemble: pd.DataFrame, complex_id: str) -> str | None:
    hits = ensemble[ensemble["complex_id"].astype(str).str.upper() == complex_id.upper()]
    if hits.empty:
        return None
    return hits.iloc[0].get("state_stratum")


def build_native_interaction_atlas(
    complexes: pd.DataFrame,
    residue_map: pd.DataFrame,
    inputs: dict[str, pd.DataFrame],
    config: dict[str, Any],
    paths: dict[str, Path],
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    benchmark = inputs["benchmark"].set_index("complex_id", drop=False)
    long_rows: list[dict[str, Any]] = []
    fp_rows: list[dict[str, Any]] = []
    native_complexes = complexes[complexes["is_native_complex"].fillna(False).astype(bool) & complexes["complex_build_status"].eq("ready")]
    for _, complex_row in native_complexes.iterrows():
        receptor_map = residue_map[residue_map["receptor_id"].astype(str).str.lower() == str(complex_row["receptor_id"]).lower()]
        interactions, meta = compute_interactions(complex_row["protein_file"], complex_row["ligand_file"], receptor_map, config)
        source = benchmark.loc[complex_row["complex_id"]] if complex_row["complex_id"] in benchmark.index else pd.Series(dtype=object)
        for idx, interaction in interactions.iterrows():
            long_rows.append(
                {
                    "native_interaction_id": f"nativeifp__{complex_row['complex_id']}__{idx}",
                    "complex_id": complex_row["complex_id"],
                    "receptor_id": complex_row["receptor_id"],
                    "ligand_id": complex_row["ligand_id"],
                    "pdb_id": complex_row["pdb_id"],
                    "ligand_comp_id": source.get("ligand_comp_id"),
                    "receptor_state": _state_for_complex(inputs["ensemble"], str(complex_row["complex_id"])),
                    "mutation_flag": source.get("mutation_flag"),
                    "mutation_list": source.get("mutation_list"),
                    "ligand_class_if_known": source.get("ligand_class"),
                    **interaction.to_dict(),
                }
            )
        bitstring, sparse_json, bits = fingerprint_from_interactions(interactions)
        fp_rows.append(
            {
                "complex_id": complex_row["complex_id"],
                "receptor_id": complex_row["receptor_id"],
                "ligand_id": complex_row["ligand_id"],
                "pdb_id": complex_row["pdb_id"],
                "receptor_state": _state_for_complex(inputs["ensemble"], str(complex_row["complex_id"])),
                "ligand_class_if_known": source.get("ligand_class"),
                "fingerprint_bitstring": bitstring,
                "fingerprint_sparse_json": sparse_json,
                "num_interactions": len(bits),
                "num_key_interactions": 0,
                "interaction_engine": meta["interaction_engine"],
                "interaction_engine_version": meta["interaction_engine_version"],
                "interaction_config_hash": meta["interaction_config_hash"],
                "warnings_json": json.dumps(meta["warnings"]),
            }
        )
    long = pd.DataFrame(long_rows, columns=NATIVE_INTERACTIONS_LONG_COLUMNS)
    fps = pd.DataFrame(fp_rows, columns=NATIVE_FINGERPRINT_COLUMNS)
    key = build_key_interaction_map(long, config, paths)
    key_bits = set(key["key_interaction_id"]) if not key.empty else set()
    if not fps.empty:
        fps["num_key_interactions"] = fps["fingerprint_sparse_json"].apply(lambda value: len(set(json.loads(value or "[]")) & key_bits))
    write_table(paths["processed"] / "native_interactions_long.parquet", long)
    write_table(paths["processed"] / "native_interactions_long.csv", long)
    write_table(paths["processed"] / "native_interaction_fingerprints.parquet", fps)
    write_table(paths["processed"] / "native_interaction_fingerprints.csv", fps)
    write_table(paths["processed"] / "key_egfr_interactions.parquet", key)
    write_table(paths["processed"] / "key_egfr_interactions.csv", key)
    return long, fps, key


def build_key_interaction_map(native_long: pd.DataFrame, config: dict[str, Any], paths: dict[str, Path] | None = None) -> pd.DataFrame:
    if native_long.empty:
        return pd.DataFrame(columns=KEY_INTERACTION_COLUMNS)
    cfg = config.get("key_interactions", {})
    freq_threshold = float(cfg.get("key_interaction_threshold_native_frequency", 0.25))
    min_count = int(cfg.get("minimum_native_count", 2))
    hydrophobic_threshold = float(cfg.get("hydrophobic_frequency_threshold", 0.5))
    total = int(native_long["complex_id"].nunique())
    rows: list[dict[str, Any]] = []
    native_long = native_long.copy()
    native_long["_key"] = native_long.apply(normalize_interaction_key, axis=1)
    manual = {(item.get("residue_role"), item.get("interaction_type")): float(item.get("weight", 1.0)) for item in cfg.get("key_interaction_manual_overrides", [])}
    for key_id, group in native_long.groupby("_key"):
        count = int(group["complex_id"].nunique())
        freq = count / total if total else 0.0
        first = group.iloc[0]
        manual_flag = (first["residue_role"], first["interaction_type"]) in manual
        threshold = hydrophobic_threshold if first["interaction_type"] == "Hydrophobic" else freq_threshold
        selected = (count >= min_count and freq >= threshold) or (manual_flag and count >= min_count)
        if not selected:
            continue
        rows.append(
            {
                "key_interaction_id": key_id,
                "receptor_state_scope": "all",
                "binding_mode_scope": "native_consensus",
                "residue_name": first["residue_name"],
                "uniprot_residue_number": first["uniprot_residue_number"],
                "klifs_position": first["klifs_position"],
                "residue_role": first["residue_role"],
                "interaction_type": first["interaction_type"],
                "native_frequency": round(freq, 3),
                "native_count": count,
                "native_total": total,
                "selection_reason": "manual_override" if manual_flag else "native_frequency",
                "weight": manual.get((first["residue_role"], first["interaction_type"]), 1.0),
                "manual_override_flag": bool(manual_flag),
                "evidence_complex_ids_json": json.dumps(sorted(group["complex_id"].dropna().astype(str).unique())),
                "warnings_json": json.dumps([]),
            }
        )
    return pd.DataFrame(rows, columns=KEY_INTERACTION_COLUMNS)
