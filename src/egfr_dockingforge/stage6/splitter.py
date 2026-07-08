from __future__ import annotations

import hashlib
from typing import Iterable

import pandas as pd
from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold

from egfr_dockingforge.common.io import write_table
from egfr_dockingforge.stage6 import schemas


def build_ranking_groups(features: pd.DataFrame, labels: pd.DataFrame, paths: dict) -> pd.DataFrame:
    merged = features[["pose_id", "docking_task_id", "ligand_id", "target_receptor_id"]].merge(
        labels[["pose_id", "rank_relevance_label"]], on="pose_id", how="left"
    )
    rows = []
    for task, group in merged.groupby("docking_task_id", sort=True):
        unique_rel = group["rank_relevance_label"].nunique(dropna=True)
        num_relevant = int(group["rank_relevance_label"].fillna(0).gt(0).sum())
        rows.append(
            {
                "group_id": task,
                "docking_task_id": task,
                "ligand_id": group["ligand_id"].iloc[0],
                "target_receptor_id": group["target_receptor_id"].iloc[0],
                "num_poses": int(len(group)),
                "num_positive_or_relevant_poses": num_relevant,
                "max_relevance_label": int(group["rank_relevance_label"].fillna(0).max()),
                "group_usable_for_ranking": bool(len(group) > 1 and unique_rel > 1),
                "exclusion_reason": "" if len(group) > 1 and unique_rel > 1 else "no_relevance_variation_or_single_pose",
            }
        )
    groups = pd.DataFrame(rows, columns=schemas.RANKING_GROUP_COLUMNS)
    write_table(paths["processed"] / "ranking_groups.parquet", groups)
    write_table(paths["processed"] / "ranking_groups.csv", groups)
    return groups


def _stable_partition(values: Iterable[str], seed: int, valid_fraction: float, test_fraction: float) -> dict[str, str]:
    mapping = {}
    for value in sorted(set(str(x) for x in values)):
        digest = hashlib.sha256(f"{seed}:{value}".encode("utf-8")).hexdigest()
        bucket = int(digest[:8], 16) / 0xFFFFFFFF
        if bucket < test_fraction:
            split = "test"
        elif bucket < test_fraction + valid_fraction:
            split = "valid"
        else:
            split = "train"
        mapping[value] = split
    if len(set(mapping.values())) < min(3, len(mapping)):
        ordered = sorted(mapping)
        for idx, value in enumerate(ordered):
            mapping[value] = ["train", "valid", "test"][idx % 3]
    return mapping


def _scaffold_ids(ligands: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for row in ligands.to_dict("records"):
        mol = None
        for key in ("prepared_ligand_file", "native_ligand_file"):
            path = row.get(key)
            if isinstance(path, str) and path:
                mol = Chem.MolFromPDBFile(path, sanitize=True, removeHs=True)
                if mol is not None:
                    break
        if mol is None:
            raise RuntimeError(f"Cannot compute Bemis-Murcko scaffold for {row.get('ligand_id')}")
        scaffold = MurckoScaffold.MurckoScaffoldSmiles(mol=mol) or f"acyclic_{row['ligand_id']}"
        rows.append({"ligand_id": row["ligand_id"], "scaffold_id": scaffold})
    return pd.DataFrame(rows)


def build_splits(features: pd.DataFrame, inputs: dict[str, pd.DataFrame], config: dict, paths: dict) -> pd.DataFrame:
    seed = int(config["splits"]["seed"])
    valid_fraction = float(config["splits"]["valid_fraction"])
    test_fraction = float(config["splits"]["test_fraction"])
    scaffolds = _scaffold_ids(inputs["ligand_prep"])
    base = features[["pose_id", "docking_task_id", "ligand_id", "target_receptor_id", "receptor_state"]].merge(scaffolds, on="ligand_id", how="left")
    base["group_id"] = base["docking_task_id"]
    rows = []
    split_specs = [
        ("ligand_holdout", "ligand_id", "all poses for held-out ligands stay together"),
        ("scaffold_holdout", "scaffold_id", "all poses for held-out Bemis-Murcko scaffolds stay together"),
        ("receptor_holdout", "target_receptor_id", "all poses for held-out receptors stay together"),
        ("state_holdout_optional", "receptor_state", "all poses for held-out receptor-state strata stay together"),
    ]
    for split_name, key, reason in split_specs:
        assignment = _stable_partition(base[key].fillna("unknown"), seed + len(split_name), valid_fraction, test_fraction)
        for row in base.to_dict("records"):
            rows.append(
                {
                    "pose_id": row["pose_id"],
                    "group_id": row["group_id"],
                    "ligand_id": row["ligand_id"],
                    "scaffold_id": row["scaffold_id"],
                    "target_receptor_id": row["target_receptor_id"],
                    "receptor_state": row["receptor_state"],
                    "split_name": split_name,
                    "split_fold": 0,
                    "train_valid_test": assignment[str(row[key] if pd.notna(row[key]) else "unknown")],
                    "split_reason": reason,
                }
            )
    splits = pd.DataFrame(rows, columns=schemas.SPLIT_COLUMNS)
    write_table(paths["processed"] / "model_splits.parquet", splits)
    write_table(paths["processed"] / "model_splits.csv", splits)
    return splits
