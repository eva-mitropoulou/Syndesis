from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pandas as pd

from egfr_dockingforge.common.io import write_table
from egfr_dockingforge.stage5.prolif_engine import compute_interactions, fingerprint_from_interactions
from egfr_dockingforge.stage5.residue_mapping import normalize_interaction_key
from egfr_dockingforge.stage5.schemas import DOCKED_FINGERPRINT_COLUMNS, DOCKED_INTERACTIONS_LONG_COLUMNS


def _write_progress(long_rows: list[dict[str, Any]], fp_rows: list[dict[str, Any]], paths: dict[str, Path]) -> None:
    long = pd.DataFrame(long_rows, columns=DOCKED_INTERACTIONS_LONG_COLUMNS)
    fps = pd.DataFrame(fp_rows, columns=DOCKED_FINGERPRINT_COLUMNS)
    write_table(paths["processed"] / "docked_pose_interactions_long.parquet", long)
    write_table(paths["processed"] / "docked_pose_interactions_long.csv", long)
    write_table(paths["processed"] / "docked_pose_fingerprints.parquet", fps)
    write_table(paths["processed"] / "docked_pose_fingerprints.csv", fps)


def compute_docked_pose_interactions(
    complexes: pd.DataFrame,
    residue_map: pd.DataFrame,
    inputs: dict[str, pd.DataFrame],
    key: pd.DataFrame,
    config: dict[str, Any],
    paths: dict[str, Path],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    scores = inputs["pose_scores"].rename(columns={"original_pose_rank": "pose_rank", "original_docking_score": "docking_score"})
    score_cols = [
        "pose_id", "docking_task_id", "ligand_id", "target_receptor_id", "task_type",
        "docking_engine", "pose_rank", "docking_score", "cnnscore", "cnnaffinity",
        "rmsd_symmetry_corrected", "stage3_pose_label", "sanity_status", "receptor_state",
    ]
    meta = scores[score_cols].set_index("pose_id", drop=False)
    key_ids = set(key["key_interaction_id"]) if not key.empty else set()
    long_rows: list[dict[str, Any]] = []
    fp_rows: list[dict[str, Any]] = []
    pose_complexes = complexes[~complexes["is_native_complex"].fillna(False).astype(bool) & complexes["complex_build_status"].eq("ready")]
    for _, complex_row in pose_complexes.iterrows():
        pose_id = complex_row["pose_id"]
        pose = meta.loc[pose_id]
        receptor_map = residue_map[residue_map["receptor_id"].astype(str).str.lower() == str(complex_row["receptor_id"]).lower()]
        interactions, engine_meta = compute_interactions(complex_row["protein_file"], complex_row["ligand_file"], receptor_map, config)
        for idx, interaction in interactions.iterrows():
            key_id = normalize_interaction_key(interaction)
            long_rows.append(
                {
                    "pose_interaction_id": f"poseifp__{pose_id}__{idx}",
                    "pose_id": pose_id,
                    "docking_task_id": pose["docking_task_id"],
                    "ligand_id": pose["ligand_id"],
                    "target_receptor_id": pose["target_receptor_id"],
                    "native_receptor_id": str(pose["ligand_id"]).rsplit("_", 1)[0],
                    "task_type": pose["task_type"],
                    "docking_engine": pose["docking_engine"],
                    "pose_rank": pose["pose_rank"],
                    "docking_score": pose["docking_score"],
                    "gnina_cnnscore": pose["cnnscore"],
                    "gnina_cnnaffinity": pose["cnnaffinity"],
                    "rmsd_symmetry_corrected": pose["rmsd_symmetry_corrected"],
                    "stage3_pose_label": pose["stage3_pose_label"],
                    "sanity_status": pose["sanity_status"],
                    "receptor_state": pose["receptor_state"],
                    **interaction.to_dict(),
                    "key_interaction_flag": key_id in key_ids,
                }
            )
        bitstring, sparse_json, bits = fingerprint_from_interactions(interactions)
        fp_rows.append(
            {
                "pose_id": pose_id,
                "docking_task_id": pose["docking_task_id"],
                "ligand_id": pose["ligand_id"],
                "target_receptor_id": pose["target_receptor_id"],
                "task_type": pose["task_type"],
                "docking_engine": pose["docking_engine"],
                "pose_rank": pose["pose_rank"],
                "fingerprint_bitstring": bitstring,
                "fingerprint_sparse_json": sparse_json,
                "num_interactions": len(bits),
                "num_key_interactions": len(bits & key_ids),
                "interaction_engine": engine_meta["interaction_engine"],
                "interaction_engine_version": engine_meta["interaction_engine_version"],
                "interaction_config_hash": engine_meta["interaction_config_hash"],
                "warnings_json": json.dumps(engine_meta["warnings"]),
            }
        )
        if len(fp_rows) % 25 == 0:
            _write_progress(long_rows, fp_rows, paths)
    _write_progress(long_rows, fp_rows, paths)
    long = pd.DataFrame(long_rows, columns=DOCKED_INTERACTIONS_LONG_COLUMNS)
    fps = pd.DataFrame(fp_rows, columns=DOCKED_FINGERPRINT_COLUMNS)
    return long, fps
