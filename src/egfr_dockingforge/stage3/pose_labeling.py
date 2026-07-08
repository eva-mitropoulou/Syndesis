from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pandas as pd

from egfr_dockingforge.common.io import write_table
from egfr_dockingforge.stage3.schemas import LABEL_COLUMNS, TASK_METRIC_COLUMNS


def label_one(rmsd: float | None, sanity_status: str, strict: float, relaxed: float) -> dict[str, Any]:
    sanity_pass = sanity_status == "pass"
    invalid = not sanity_pass
    strict_flag = rmsd is not None and rmsd <= strict and sanity_pass
    relaxed_flag = rmsd is not None and rmsd <= relaxed and sanity_pass
    if invalid:
        label = "invalid_pose"
    elif strict_flag:
        label = "strict_native_like"
    elif relaxed_flag:
        label = "relaxed_native_like"
    else:
        label = "pending_review"
    return {"strict": strict_flag, "relaxed": relaxed_flag, "invalid": invalid, "label": label}


def build_labels(rmsd: pd.DataFrame, sanity: pd.DataFrame, config: dict[str, Any], paths: dict[str, Path]) -> pd.DataFrame:
    if rmsd.empty:
        frame = pd.DataFrame(columns=LABEL_COLUMNS)
    else:
        sanity_lookup = sanity.set_index("pose_id")
        rows = []
        for _, row in rmsd.iterrows():
            s = sanity_lookup.loc[row["pose_id"]]
            flags = label_one(row["rmsd_symmetry_corrected"], s["sanity_status"], config["rmsd"]["strict_threshold_angstrom"], config["rmsd"]["relaxed_threshold_angstrom"])
            rows.append({
                "pose_id": row["pose_id"],
                "docking_task_id": row["docking_task_id"],
                "ligand_id": row["ligand_id"],
                "target_receptor_id": row["target_receptor_id"],
                "task_type": row["task_type"],
                "docking_engine": row["docking_engine"],
                "pose_rank": row["pose_rank"],
                "docking_score": row["docking_score"],
                "rmsd_symmetry_corrected": row["rmsd_symmetry_corrected"],
                "sanity_status": s["sanity_status"],
                "strict_native_like_flag": flags["strict"],
                "relaxed_native_like_flag": flags["relaxed"],
                "invalid_pose_flag": flags["invalid"],
                "sampled_not_ranked_flag": False,
                "ranking_failure_flag": False,
                "sampling_failure_flag": False,
                "stage3_pose_label": flags["label"],
                "interaction_recovery_status": "pending_stage5",
                "final_pose_label_status": "pending_stage5",
                "label_warnings_json": json.dumps([]),
            })
        frame = pd.DataFrame(rows, columns=LABEL_COLUMNS)
        top_n = int(config["rmsd"].get("top_n_sampling", 20))
        for task_id, group in frame.groupby("docking_task_id"):
            top = group.sort_values("pose_rank")
            top1_idx = top.index[0]
            topn = top[top["pose_rank"] <= top_n]
            relaxed_exists = bool(topn["relaxed_native_like_flag"].any())
            strict_exists = bool(topn["strict_native_like_flag"].any())
            top1_relaxed = bool(frame.loc[top1_idx, "relaxed_native_like_flag"])
            top1_strict = bool(frame.loc[top1_idx, "strict_native_like_flag"])
            if relaxed_exists and not top1_relaxed:
                idxs = topn.index[topn["relaxed_native_like_flag"]]
                frame.loc[idxs, "sampled_not_ranked_flag"] = True
                wrong = top.index[~top["relaxed_native_like_flag"] & ~top["invalid_pose_flag"]]
                frame.loc[wrong, "ranking_failure_flag"] = True
                frame.loc[wrong, "stage3_pose_label"] = "ranking_failure"
            if not relaxed_exists:
                idxs = top.index[~top["invalid_pose_flag"]]
                frame.loc[idxs, "sampling_failure_flag"] = True
                frame.loc[idxs, "stage3_pose_label"] = "sampling_failure"
            if strict_exists and not top1_strict and top1_relaxed:
                frame.loc[top1_idx, "stage3_pose_label"] = "relaxed_native_like"
    write_table(paths["processed"] / "stage3_pose_labels.parquet", frame)
    write_table(paths["processed"] / "stage3_pose_labels.csv", frame)
    return frame


def task_metrics(tasks: pd.DataFrame, labels: pd.DataFrame, runs: pd.DataFrame, config: dict[str, Any], paths: dict[str, Path]) -> pd.DataFrame:
    top_n = int(config["rmsd"].get("top_n_sampling", 20))
    rows = []
    label_groups = {k: v for k, v in labels.groupby("docking_task_id")} if not labels.empty else {}
    run_lookup = runs.groupby("docking_task_id")["status"].first().to_dict() if not runs.empty else {}
    run_success = {key: value == "success" for key, value in run_lookup.items()}
    for _, task in tasks.iterrows():
        group = label_groups.get(task["docking_task_id"], pd.DataFrame())
        failed = not run_success.get(task["docking_task_id"], False)
        sorted_group = group.sort_values("pose_rank") if not group.empty else group
        top1 = sorted_group.iloc[0] if not sorted_group.empty else None
        def best(limit: int) -> float | None:
            if sorted_group.empty:
                return None
            vals = sorted_group[sorted_group["pose_rank"] <= limit]["rmsd_symmetry_corrected"].dropna()
            return float(vals.min()) if not vals.empty else None
        def first_rank(column: str) -> int | None:
            if sorted_group.empty:
                return None
            hits = sorted_group[sorted_group[column] == True]
            return int(hits["pose_rank"].min()) if not hits.empty else None
        sanity_pass = None
        if not sorted_group.empty:
            sanity_pass = float((sorted_group["sanity_status"] == "pass").mean())
        relaxed_success = bool((sorted_group[sorted_group["pose_rank"] <= top_n]["relaxed_native_like_flag"] == True).any()) if not sorted_group.empty else False
        strict_success = bool((sorted_group[sorted_group["pose_rank"] <= top_n]["strict_native_like_flag"] == True).any()) if not sorted_group.empty else False
        top1_relaxed = bool(top1 is not None and top1["relaxed_native_like_flag"])
        top1_strict = bool(top1 is not None and top1["strict_native_like_flag"])
        if failed:
            failure_category = run_lookup.get(task["docking_task_id"], "failed_run")
        elif relaxed_success and not top1_relaxed:
            failure_category = "ranking_failure"
        elif not relaxed_success:
            failure_category = "sampling_failure"
        else:
            failure_category = "success"
        rows.append({
            "docking_task_id": task["docking_task_id"],
            "task_type": task["task_type"],
            "ligand_id": task["ligand_id"],
            "native_receptor_id": task["native_receptor_id"],
            "target_receptor_id": task["target_receptor_id"],
            "state_match_flag": task["state_match_flag"],
            "docking_engine": task["docking_engine"],
            "num_poses_generated": int(len(group)),
            "top1_rmsd": None if top1 is None else top1["rmsd_symmetry_corrected"],
            "top1_score": None if top1 is None else top1["docking_score"],
            "top5_best_rmsd": best(5),
            "top10_best_rmsd": best(10),
            "top20_best_rmsd": best(20),
            "best_rmsd_any_pose": best(9999),
            "rank_of_first_strict_native_like": first_rank("strict_native_like_flag"),
            "rank_of_first_relaxed_native_like": first_rank("relaxed_native_like_flag"),
            "sampling_success_strict_top20": strict_success,
            "sampling_success_relaxed_top20": relaxed_success,
            "ranking_success_strict_top1": top1_strict,
            "ranking_success_relaxed_top1": top1_relaxed,
            "physical_sanity_pass_rate": sanity_pass,
            "failed_run_flag": failed,
            "failure_category": failure_category,
            "warnings_json": json.dumps(["No docked poses available."] if failed else []),
        })
    frame = pd.DataFrame(rows, columns=TASK_METRIC_COLUMNS)
    write_table(paths["processed"] / "docking_task_metrics.parquet", frame)
    write_table(paths["processed"] / "docking_task_metrics.csv", frame)
    return frame
