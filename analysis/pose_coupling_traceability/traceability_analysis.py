#!/usr/bin/env python3
"""Audit whether independently maximized late-fusion terms share a pose.

This is a deterministic, post hoc analysis of frozen receptor-level outputs.  It
does not dock, rescore, fingerprint, bootstrap, or permute any molecule.
"""
from __future__ import annotations

import argparse
import json
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
OUT = Path(__file__).resolve().parent
RTOL = 1e-10
ATOL = 1e-12


def packaged_path(public_relative: str, working_relative: str) -> Path:
    """Select the compact public-package input when it is available."""
    public = ROOT / public_relative
    return public if public.exists() else ROOT / working_relative


@dataclass(frozen=True)
class TargetSpec:
    target: str
    master: Path
    receptors: tuple[str, ...]
    expected_n: int
    expected_actives: int
    expected_top_counts: tuple[int, int, int | None]


SPECS = (
    TargetSpec(
        "EGFR",
        packaged_path(
            "data/benchmark/egfr_pose_scores.parquet",
            "results_showcase/submission_robustness/corrected_enrichment/egfr_master.parquet",
        ),
        ("1m17_a_aq4_999", "1xkk_a_fmm_91", "4hjo_a_aq4_1001", "5cav_a_4zq_1101"),
        35_552,
        542,
        (65, 89, 86),
    ),
    TargetSpec(
        "CDK2",
        packaged_path(
            "data/benchmark/cdk2_pose_scores.parquet",
            "results_showcase/submission_robustness/cdk2_four_receptor_four_prior/cdk2_master.parquet",
        ),
        ("1fin_a_atp", "2a4l_a_rrc", "1aq1_a_stu", "1pxn_a_ck6"),
        28_296,
        474,
        (55, 71, None),
    ),
)


def read_parquet(path: Path) -> pd.DataFrame:
    """Use the declared PyArrow stack when available, with a local fallback."""
    try:
        return pd.read_parquet(path)
    except ImportError:
        return pd.read_parquet(path, engine="fastparquet")


def parse_bits(value: object) -> set[str]:
    return set(json.loads(value)) if isinstance(value, str) and value else set()


def top_membership(scores: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Stable descending ranks and the manuscript top-1% membership mask."""
    order = np.argsort(-scores, kind="mergesort")
    ranks = np.empty(len(scores), dtype=int)
    ranks[order] = np.arange(1, len(scores) + 1)
    top = np.zeros(len(scores), dtype=bool)
    top[order[: max(1, round(0.01 * len(scores)))]] = True
    return ranks, top


def maxima(values: np.ndarray, receptors: tuple[str, ...]) -> tuple[list[str], str]:
    maximum = np.nanmax(values)
    members = [receptor for receptor, value in zip(receptors, values) if np.isfinite(value) and np.isclose(value, maximum, rtol=RTOL, atol=ATOL)]
    if not members:
        raise AssertionError("A ligand had no finite receptor-level maximum")
    return members, members[0]


def summarise_group(frame: pd.DataFrame, target: str, scope: str, group: str, mask: pd.Series) -> dict[str, object]:
    subset = frame.loc[mask]
    n = len(subset)
    nonreal = int((~subset["late_fusion_pose_realizable"]).sum())
    return {
        "target": target,
        "coverage_scope": scope,
        "group": group,
        "n_molecules": n,
        "pose_realizable_count": n - nonreal,
        "pose_nonrealizable_count": nonreal,
        "pose_nonrealizable_percent": 100 * nonreal / n if n else np.nan,
    }


def gap_summary(frame: pd.DataFrame, target: str, scope: str, group: str, mask: pd.Series) -> dict[str, object]:
    values = frame.loc[mask, "fusion_gap"].to_numpy(float)
    if len(values) == 0:
        return {"target": target, "coverage_scope": scope, "group": group, "n_molecules": 0}
    return {
        "target": target,
        "coverage_scope": scope,
        "group": group,
        "n_molecules": len(values),
        "median": float(np.median(values)),
        "iqr_q1": float(np.percentile(values, 25)),
        "iqr_q3": float(np.percentile(values, 75)),
        "p95": float(np.percentile(values, 95)),
        "maximum": float(np.max(values)),
        "percent_gt_zero": float(100 * np.mean(values > ATOL)),
    }


def prepare_target(spec: TargetSpec) -> tuple[pd.DataFrame, pd.DataFrame]:
    raw = read_parquet(spec.master).copy()
    required = {"lit_pcba_id", "label", "target_receptor_id", "cnnscore", "key_interaction_recall_consensus", "fingerprint_sparse_json", "status"}
    missing = required.difference(raw.columns)
    if missing:
        raise RuntimeError(f"{spec.target}: master is missing required columns: {sorted(missing)}")
    raw["ligand_id"] = raw["lit_pcba_id"].astype(str)
    raw = raw[raw["target_receptor_id"].isin(spec.receptors)].copy()
    if raw.duplicated(["ligand_id", "target_receptor_id"]).any():
        raise RuntimeError(f"{spec.target}: duplicate ligand-receptor rows")
    counts = raw.groupby("target_receptor_id")["ligand_id"].nunique().to_dict()
    if any(counts.get(receptor) != spec.expected_n for receptor in spec.receptors):
        raise RuntimeError(f"{spec.target}: incomplete receptor table: {counts}")
    if raw["label"].isna().any() or raw.groupby("ligand_id")["label"].nunique().gt(1).any():
        raise RuntimeError(f"{spec.target}: labels are missing or inconsistent")

    # The frozen manuscript analyses obtain the ligand index through pivot(),
    # which is lexical by ligand identifier.  Retain that stable order here.
    ligand_ids = np.array(sorted(raw["ligand_id"].unique()), dtype=object)
    labels = raw.drop_duplicates("ligand_id").set_index("ligand_id")["label"].reindex(ligand_ids).to_numpy(int)
    if len(ligand_ids) != spec.expected_n or int(labels.sum()) != spec.expected_actives:
        raise RuntimeError(f"{spec.target}: unexpected molecule/active counts")

    receptor_index = pd.MultiIndex.from_product([ligand_ids, spec.receptors], names=["ligand_id", "target_receptor_id"])
    indexed = raw.set_index(["ligand_id", "target_receptor_id"]).reindex(receptor_index)
    cnn_matrix = indexed["cnnscore"].to_numpy(float).reshape(len(ligand_ids), len(spec.receptors))
    recall_matrix = indexed["key_interaction_recall_consensus"].to_numpy(float).reshape(len(ligand_ids), len(spec.receptors))
    coupled_matrix = cnn_matrix * (1.0 + recall_matrix)
    rows: list[dict[str, object]] = []
    for index, ligand_id in enumerate(ligand_ids):
        cnn = cnn_matrix[index]
        recall = recall_matrix[index]
        available = np.isfinite(cnn) & np.isfinite(recall)
        if not available.any():
            raise RuntimeError(f"{spec.target}: {ligand_id} has no available receptor-specific pose")
        coupled = coupled_matrix[index]
        cnn_set, cnn_first = maxima(cnn, spec.receptors)
        recall_set, recall_first = maxima(recall, spec.receptors)
        coupled_set, coupled_first = maxima(coupled, spec.receptors)
        realizable = bool(set(cnn_set).intersection(recall_set))
        gnina_score = float(np.nanmax(cnn))
        coupled_score = float(np.nanmax(coupled))
        late_score = float(gnina_score * (1.0 + np.nanmax(recall)))
        gap = late_score - coupled_score
        if gap < -ATOL:
            raise AssertionError(f"{spec.target}: negative fusion gap for {ligand_id}: {gap}")
        rows.append({
            "ligand_id": ligand_id,
            "target": spec.target,
            "activity_label": int(labels[index]),
            "activity": "active" if labels[index] else "decoy",
            "stable_ligand_order": index + 1,
            "gnina_score": gnina_score,
            "coupled_score": coupled_score,
            "late_fusion_score": late_score,
            "max_cnn_receptor_set": ";".join(cnn_set),
            "max_recall_receptor_set": ";".join(recall_set),
            "max_coupled_receptor_set": ";".join(coupled_set),
            "max_cnn_receptor_stable": cnn_first,
            "max_recall_receptor_stable": recall_first,
            "max_coupled_receptor_stable": coupled_first,
            "late_fusion_pose_realizable": realizable,
            "fusion_gap": max(0.0, gap),
            "receptor_coverage_count": int(available.sum()),
        })
    per_ligand = pd.DataFrame(rows)
    for arm, score_column in (("gnina", "gnina_score"), ("coupled", "coupled_score"), ("late_fusion", "late_fusion_score")):
        rank, top = top_membership(per_ligand[score_column].to_numpy(float))
        per_ligand[f"{arm}_rank"] = rank
        per_ligand[f"{arm}_top1"] = top
    per_ligand["rescued_active"] = per_ligand["activity_label"].eq(1) & per_ligand["coupled_top1"] & ~per_ligand["gnina_top1"]
    per_ligand["lost_active"] = per_ligand["activity_label"].eq(1) & per_ligand["gnina_top1"] & ~per_ligand["coupled_top1"]
    per_ligand["shared_active"] = per_ligand["activity_label"].eq(1) & per_ligand["gnina_top1"] & per_ligand["coupled_top1"]

    top_counts = tuple(int(per_ligand.loc[per_ligand[f"{arm}_top1"], "activity_label"].sum()) for arm in ("gnina", "coupled", "late_fusion"))
    for observed, expected, arm in zip(top_counts, spec.expected_top_counts, ("GNINA", "coupled", "late fusion")):
        if expected is not None and observed != expected:
            raise AssertionError(f"{spec.target}: {arm} top-1% active count {observed} != published {expected}")
    if not np.allclose(per_ligand["late_fusion_score"], per_ligand["gnina_score"] * (1.0 + np.nanmax(recall_matrix, axis=1)), rtol=RTOL, atol=ATOL):
        raise AssertionError(f"{spec.target}: late-fusion formula check failed")
    if "pose_file" not in indexed.columns:
        indexed["pose_file"] = ""
    receptor_rows = indexed.reset_index()[["ligand_id", "target_receptor_id", "pose_file", "fingerprint_sparse_json"]].rename(columns={"target_receptor_id": "receptor"})
    receptor_rows["target"] = spec.target
    receptor_rows["activity_label"] = np.repeat(labels, len(spec.receptors))
    receptor_rows["available"] = np.isfinite(cnn_matrix).ravel() & np.isfinite(recall_matrix).ravel()
    receptor_rows["cnnscore"] = cnn_matrix.ravel()
    receptor_rows["native_interaction_recall"] = recall_matrix.ravel()
    receptor_rows["coupled_score"] = coupled_matrix.ravel()
    receptor_rows["pose_file"] = receptor_rows["pose_file"].fillna("")
    receptor_rows["fingerprint_sparse_json"] = receptor_rows["fingerprint_sparse_json"].fillna("")
    return per_ligand, receptor_rows


def group_outputs(per_ligand: pd.DataFrame, target: str) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    subset_rows: list[dict[str, object]] = []
    gap_rows: list[dict[str, object]] = []
    activity_rows: list[dict[str, object]] = []
    scopes = [("all_available", pd.Series(True, index=per_ligand.index))]
    if target == "CDK2":
        scopes.append(("complete_four_receptor", per_ligand["receptor_coverage_count"].eq(4)))
    for scope, coverage in scopes:
        local_top: dict[str, pd.Series] = {}
        for arm, score_column in (("gnina", "gnina_score"), ("coupled", "coupled_score"), ("late_fusion", "late_fusion_score")):
            if scope == "all_available":
                local_top[arm] = per_ligand[f"{arm}_top1"]
            else:
                values = per_ligand.loc[coverage, score_column].to_numpy(float)
                _, membership = top_membership(values)
                mask = pd.Series(False, index=per_ligand.index)
                mask.loc[per_ligand.index[coverage]] = membership
                local_top[arm] = mask
        groups = {
            "all_molecules": coverage,
            "all_actives": coverage & per_ligand["activity_label"].eq(1),
            "all_decoys": coverage & per_ligand["activity_label"].eq(0),
            "gnina_top1": coverage & local_top["gnina"],
            "coupled_top1": coverage & local_top["coupled"],
            "late_fusion_top1": coverage & local_top["late_fusion"],
            "actives_in_coupled_top1": coverage & per_ligand["activity_label"].eq(1) & local_top["coupled"],
            "actives_in_late_fusion_top1": coverage & per_ligand["activity_label"].eq(1) & local_top["late_fusion"],
        }
        subset_rows.extend(summarise_group(per_ligand, target, scope, name, mask) for name, mask in groups.items())
        gap_groups = {
            "all_molecules": coverage,
            "actives": coverage & per_ligand["activity_label"].eq(1),
            "coupled_top1": coverage & local_top["coupled"],
            "rescued_actives": coverage & per_ligand["rescued_active"],
            "pose_nonrealizable_only": coverage & ~per_ligand["late_fusion_pose_realizable"],
        }
        gap_rows.extend(gap_summary(per_ligand, target, scope, name, mask) for name, mask in gap_groups.items())

    for name, mask in {
        "rescued_actives": per_ligand["rescued_active"],
        "lost_actives": per_ligand["lost_active"],
        "shared_actives": per_ligand["shared_active"],
    }.items():
        subset = per_ligand.loc[mask]
        activity_rows.append({
            "target": target,
            "group": name,
            "n_actives": len(subset),
            "pose_nonrealizable_count": int((~subset["late_fusion_pose_realizable"]).sum()),
            "pose_nonrealizable_percent": 100 * float((~subset["late_fusion_pose_realizable"]).mean()) if len(subset) else np.nan,
            "median_gnina_rank": float(subset["gnina_rank"].median()) if len(subset) else np.nan,
            "median_coupled_rank": float(subset["coupled_rank"].median()) if len(subset) else np.nan,
            "median_late_fusion_rank": float(subset["late_fusion_rank"].median()) if len(subset) else np.nan,
            "median_gnina_to_coupled_rank_change": float((subset["gnina_rank"] - subset["coupled_rank"]).median()) if len(subset) else np.nan,
        })
    return pd.DataFrame(subset_rows), pd.DataFrame(gap_rows), pd.DataFrame(activity_rows)


def tie_and_transition_outputs(per_ligand: pd.DataFrame, target: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    rows: list[dict[str, object]] = []
    scopes = [("all_available", pd.Series(True, index=per_ligand.index))]
    if target == "CDK2":
        scopes.append(("complete_four_receptor", per_ligand["receptor_coverage_count"].eq(4)))
    for scope, mask in scopes:
        subset = per_ligand.loc[mask]
        record: dict[str, object] = {"target": target, "coverage_scope": scope, "n_molecules": len(subset)}
        for name, column in (("cnn", "max_cnn_receptor_set"), ("recall", "max_recall_receptor_set"), ("coupled", "max_coupled_receptor_set")):
            ties = subset[column].str.count(";").add(1).gt(1)
            record[f"unique_max_{name}_count"] = int((~ties).sum())
            record[f"tied_max_{name}_count"] = int(ties.sum())
        rows.append(record)
    rescued = per_ligand.loc[per_ligand["rescued_active"]].copy()
    transitions = (
        rescued.groupby(["target", "max_cnn_receptor_stable", "max_recall_receptor_stable", "max_coupled_receptor_stable"], as_index=False)
        .size().rename(columns={"size": "n_rescued_actives"})
        .sort_values(["target", "n_rescued_actives", "max_cnn_receptor_stable", "max_recall_receptor_stable", "max_coupled_receptor_stable"], ascending=[True, False, True, True, True], kind="mergesort")
    )
    transitions["percent_of_target_rescued_actives"] = 100 * transitions["n_rescued_actives"] / transitions.groupby("target")["n_rescued_actives"].transform("sum")
    return pd.DataFrame(rows), transitions


def egfr_native_union() -> set[str]:
    native = read_parquet(packaged_path(
        "data/benchmark/egfr_native_fingerprints.parquet",
        "data/processed/stage5/native_interaction_fingerprints.parquet",
    ))
    primary = set(SPECS[0].receptors)
    native = native[native["receptor_id"].isin(primary)]
    union: set[str] = set()
    for value in native["fingerprint_sparse_json"]:
        union.update(parse_bits(value))
    if len(union) != 62:
        raise AssertionError(f"EGFR native union has {len(union)} bits, expected 62")
    return union


def representative_case(per_ligand: pd.DataFrame, receptor_rows: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    candidates = per_ligand[
        per_ligand["activity_label"].eq(1)
        & per_ligand["coupled_top1"]
        & ~per_ligand["gnina_top1"]
        & ~per_ligand["late_fusion_pose_realizable"]
    ].copy()
    if candidates.empty:
        raise RuntimeError("No active, coupled-top-1%, GNINA-rescued, pose-non-realizable EGFR molecule exists")
    candidates["coupled_differs_from_both"] = (
        candidates["max_coupled_receptor_stable"].ne(candidates["max_cnn_receptor_stable"])
        & candidates["max_coupled_receptor_stable"].ne(candidates["max_recall_receptor_stable"])
    )
    if candidates["coupled_differs_from_both"].any():
        candidates = candidates[candidates["coupled_differs_from_both"]].copy()
    candidates["rank_improvement"] = candidates["gnina_rank"] - candidates["coupled_rank"]
    candidates = candidates.sort_values(["rank_improvement", "stable_ligand_order"], ascending=[False, True], kind="mergesort")
    selected = candidates.iloc[0]
    rows = receptor_rows[receptor_rows["ligand_id"].eq(selected.ligand_id)].copy()
    rows["cnnscore_rank_among_receptors"] = rows["cnnscore"].rank(ascending=False, method="min")
    rows["recall_rank_among_receptors"] = rows["native_interaction_recall"].rank(ascending=False, method="min")
    rows["coupled_score_rank_among_receptors"] = rows["coupled_score"].rank(ascending=False, method="min")
    rows["selected_by_gnina"] = rows["receptor"].isin(str(selected.max_cnn_receptor_set).split(";"))
    rows["selected_by_recall_maximum"] = rows["receptor"].isin(str(selected.max_recall_receptor_set).split(";"))
    rows["selected_by_coupled_score"] = rows["receptor"].isin(str(selected.max_coupled_receptor_set).split(";"))
    union = egfr_native_union()
    interaction_rows: list[dict[str, object]] = []
    for row in rows.itertuples():
        recovered = sorted(parse_bits(row.fingerprint_sparse_json).intersection(union))
        rows.loc[rows["receptor"].eq(row.receptor), "recovered_native_interaction_bit_count"] = len(recovered)
        for interaction in recovered:
            interaction_rows.append({"ligand_id": selected.ligand_id, "receptor": row.receptor, "interaction_bit": interaction})
    representative = rows.merge(selected.to_frame().T, on=["ligand_id", "target", "activity_label"], how="left", suffixes=("", "_ligand"))
    return selected.to_frame().T, representative, pd.DataFrame(interaction_rows)


def write_pymol_and_figure(representative: pd.DataFrame, interactions: pd.DataFrame) -> None:
    """Copy the three relevant poses, write a portable PML, and render a compact PNG."""
    roles = (
        ("maximum CNNscore", representative["selected_by_gnina"]),
        ("maximum recall", representative["selected_by_recall_maximum"]),
        ("coupled-score selection", representative["selected_by_coupled_score"]),
    )
    pose_dir = OUT / "representative_poses"
    pose_dir.mkdir(exist_ok=True)
    pml_lines = ["reinitialize", "bg_color white", "set ray_opaque_background, off", "set cartoon_transparency, 0.15"]
    panel_paths: list[tuple[str, Path, pd.Series]] = []
    receptor_dir = ROOT / "data/processed/stage3/docking_receptors"
    for number, (role, mask) in enumerate(roles, start=1):
        matches = representative.loc[mask]
        if len(matches) != 1:
            raise AssertionError(f"Representative {role} does not have one deterministic receptor")
        row = matches.iloc[0]
        source = Path(row.pose_file)
        if not source.exists():
            raise RuntimeError(f"Representative pose file is unavailable: {source}")
        pose_name = f"{row.receptor}_{row.ligand_id}.pdbqt"
        destination = pose_dir / pose_name
        if not destination.exists():
            shutil.copy2(source, destination)
        receptor = receptor_dir / f"{row.receptor}.pdbqt"
        if not receptor.exists():
            raise RuntimeError(f"Representative receptor file is unavailable: {receptor}")
        object_prefix = f"panel_{number}"
        pml_lines.extend([
            f"# Panel {chr(64 + number)}: {role}",
            f"load representative_poses/receptors/{row.receptor}.pdbqt, {object_prefix}_receptor",
            f"load representative_poses/{pose_name}, {object_prefix}_ligand",
            f"hide everything, {object_prefix}_receptor",
            f"show cartoon, {object_prefix}_receptor",
            f"color gray80, {object_prefix}_receptor",
            f"show sticks, {object_prefix}_ligand",
            f"util.cbag {object_prefix}_ligand",
            f"select {object_prefix}_pocket, byres ({object_prefix}_receptor within 5 of {object_prefix}_ligand)",
            f"show sticks, {object_prefix}_pocket",
            f"color slate, {object_prefix}_pocket",
            f"zoom {object_prefix}_ligand, 12",
            f"save representative_{number}_{row.receptor}.pse",
        ])
        panel_pml = OUT / f"representative_panel_{number}.pml"
        panel_pml.write_text("\n".join(pml_lines[-12:] + [f"png representative_panel_{number}.png, 1000, 800, dpi=220, ray=1"]) + "\n")
        panel_paths.append((f"Panel {chr(64 + number)} — {role}: {row.receptor}", OUT / f"representative_panel_{number}.png", row))
    (OUT / "representative_case.pml").write_text("\n".join(pml_lines) + "\n")

    pymol = shutil.which("pymol")
    if not pymol:
        return
    try:
        for number, (_, _, row) in enumerate(panel_paths, start=1):
            subprocess.run([pymol, "-cq", str(OUT / f"representative_panel_{number}.pml")], cwd=OUT, check=True, capture_output=True, text=True)
    except subprocess.CalledProcessError:
        (OUT / "rendering_status.md").write_text(
            "PyMOL was detected but could not run in the current Python environment. "
            "The portable representative_case.pml and the copied pose files are provided for reproducible rendering.\n"
        )
        return
    import matplotlib.pyplot as plt
    from matplotlib.image import imread
    figure, axes = plt.subplots(1, len(panel_paths), figsize=(5.0 * len(panel_paths), 5.8), constrained_layout=True)
    if len(panel_paths) == 1:
        axes = [axes]
    for axis, (title, image_path, row) in zip(axes, panel_paths):
        axis.imshow(imread(image_path))
        axis.set_title(title.replace("_", " "), fontsize=10, fontweight="bold")
        key = interactions[interactions["receptor"].eq(row.receptor)]["interaction_bit"].tolist()
        non_vdw = [bit for bit in key if not bit.endswith(":VdWContact")]
        shown = (non_vdw or key)[:5]
        text = (
            f"CNNscore {row.cnnscore:.5f} | recall {row.native_interaction_recall:.3f} | coupled {row.coupled_score:.5f}\n"
            f"Recovered native bits ({int(row.recovered_native_interaction_bit_count)}): " + ", ".join(shown)
        )
        axis.text(0.01, -0.05, text, transform=axis.transAxes, ha="left", va="top", fontsize=7, wrap=True)
        axis.axis("off")
    figure.savefig(OUT / "representative_case.png", dpi=300, bbox_inches="tight")
    plt.close(figure)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--render", action="store_true", help="Render the representative figure with local source pose files and PyMOL.")
    args = parser.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    all_ligands, all_subsets, all_gaps, all_activity, all_ties, all_transitions = [], [], [], [], [], []
    egfr_ligands = egfr_receptors = None
    for spec in SPECS:
        per_ligand, receptor_rows = prepare_target(spec)
        subsets, gaps, activity = group_outputs(per_ligand, spec.target)
        ties, transitions = tie_and_transition_outputs(per_ligand, spec.target)
        all_ligands.append(per_ligand)
        all_subsets.append(subsets)
        all_gaps.append(gaps)
        all_activity.append(activity)
        all_ties.append(ties)
        all_transitions.append(transitions)
        if spec.target == "EGFR":
            egfr_ligands, egfr_receptors = per_ligand, receptor_rows
    assert egfr_ligands is not None and egfr_receptors is not None
    per_ligand = pd.concat(all_ligands, ignore_index=True)
    subsets = pd.concat(all_subsets, ignore_index=True)
    gaps = pd.concat(all_gaps, ignore_index=True)
    activity = pd.concat(all_activity, ignore_index=True)
    ties = pd.concat(all_ties, ignore_index=True)
    transitions = pd.concat(all_transitions, ignore_index=True)
    selected, representative, interactions = representative_case(egfr_ligands, egfr_receptors)

    # Validation assertions requested for a reusable frozen-data analysis.
    assert not per_ligand.duplicated(["target", "ligand_id"]).any()
    assert all(per_ligand.groupby("target").size().reindex([s.target for s in SPECS]).to_numpy() == [s.expected_n for s in SPECS])
    assert not (per_ligand["fusion_gap"] < -ATOL).any()
    assert per_ligand.loc[per_ligand["target"].eq("CDK2") & per_ligand["receptor_coverage_count"].eq(3), "activity_label"].eq(0).all()
    assert bool(selected["rescued_active"].iloc[0]) and not bool(selected["late_fusion_pose_realizable"].iloc[0])

    per_ligand.to_csv(OUT / "traceability_by_ligand.csv", index=False)
    subsets.to_csv(OUT / "traceability_subsets.csv", index=False)
    gaps.to_csv(OUT / "fusion_gap_summary.csv", index=False)
    activity.to_csv(OUT / "rescued_lost_shared_actives.csv", index=False)
    per_ligand.loc[per_ligand["rescued_active"]].to_csv(OUT / "rescued_actives.csv", index=False)
    selected.to_csv(OUT / "representative_selection.csv", index=False)
    representative.drop(columns=["pose_file", "fingerprint_sparse_json"]).to_csv(OUT / "representative_case.csv", index=False)
    interactions.to_csv(OUT / "representative_case_interactions.csv", index=False)
    summary = subsets[subsets["group"].eq("all_molecules")].merge(
        gaps[gaps["group"].eq("all_molecules")][["target", "coverage_scope", "median", "iqr_q1", "iqr_q3", "p95", "maximum", "percent_gt_zero"]],
        on=["target", "coverage_scope"], how="left",
    )
    summary = summary.merge(ties.drop(columns="n_molecules"), on=["target", "coverage_scope"], how="left")
    summary.to_csv(OUT / "traceability_summary.csv", index=False)
    ties.to_csv(OUT / "maxima_tie_summary.csv", index=False)
    transitions.to_csv(OUT / "rescued_active_receptor_transitions.csv", index=False)
    if args.render:
        write_pymol_and_figure(representative, interactions)
    print(summary.to_string(index=False))
    print(activity.to_string(index=False))
    print(selected[["ligand_id", "gnina_rank", "coupled_rank", "late_fusion_rank", "fusion_gap"]].to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
