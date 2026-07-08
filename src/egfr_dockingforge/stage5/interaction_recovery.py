from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pandas as pd

from egfr_dockingforge.common.io import write_table
from egfr_dockingforge.stage5.prolif_engine import tanimoto
from egfr_dockingforge.stage5.schemas import INTERACTION_RECOVERY_COLUMNS


def _bits(value: object) -> set[str]:
    try:
        return set(json.loads(value if isinstance(value, str) and value else "[]"))
    except json.JSONDecodeError:
        return set()


def _metrics(pose_bits: set[str], ref_bits: set[str], key_bits: set[str]) -> tuple[float | None, float | None, float | None, float | None, set[str], set[str], set[str]]:
    pose_key = pose_bits & key_bits
    ref_key = ref_bits & key_bits
    recovered = pose_key & ref_key
    missing = ref_key - pose_key
    extra = pose_key - ref_key
    recall = None if not ref_key else len(recovered) / len(ref_key)
    precision = None if not pose_key else len(recovered) / len(pose_key)
    if recall is None or precision is None or recall + precision == 0:
        f1 = None
    else:
        f1 = 2 * recall * precision / (recall + precision)
    return tanimoto(pose_bits, ref_bits), recall, precision, f1, missing, extra, recovered


def _label(recall: float | None, tanimoto_value: float | None, config: dict[str, Any]) -> str:
    cfg = config.get("recovery", {})
    if recall is None or tanimoto_value is None:
        return "no_reference"
    if recall >= float(cfg.get("high_key_recall_threshold", 0.75)) and tanimoto_value >= float(cfg.get("high_tanimoto_threshold", 0.5)):
        return "high_recovery"
    if recall >= float(cfg.get("moderate_key_recall_threshold", 0.5)):
        return "moderate_recovery"
    return "poor_recovery"


def compute_interaction_recovery(
    docked_fps: pd.DataFrame,
    native_fps: pd.DataFrame,
    key: pd.DataFrame,
    inputs: dict[str, pd.DataFrame],
    config: dict[str, Any],
    paths: dict[str, Path],
) -> pd.DataFrame:
    scores = inputs["pose_scores"].rename(columns={"original_pose_rank": "pose_rank"})
    score_meta = scores.set_index("pose_id", drop=False)
    native_by_receptor = {str(row.receptor_id).lower(): row for row in native_fps.itertuples(index=False)}
    consensus_bits = set().union(*[_bits(value) for value in native_fps["fingerprint_sparse_json"]]) if not native_fps.empty else set()
    # Binding-mode recovery is scored against the CONSERVED CORE of key
    # interactions (those present in a large fraction of native complexes), not
    # the union of every native's idiosyncratic contacts. Scoring recall against
    # the union makes even a co-crystal ligand in its own receptor score ~0.2-0.4
    # (most "key" bits belong to other, chemically different ligands), which
    # makes the downstream binding-mode gate impossible to pass. The conserved
    # core (hinge, gatekeeper, catalytic Lys/Glu, DFG) is the pharmacophore we
    # actually want analogs to preserve.
    core_min_freq = float(config.get("recovery", {}).get("consensus_core_min_frequency", 0.6))
    if not key.empty and "native_frequency" in key.columns:
        core_key = key[
            (key["native_frequency"].astype(float) >= core_min_freq)
            | key["manual_override_flag"].fillna(False).astype(bool)
        ]
        key_bits = set(core_key["key_interaction_id"]) if not core_key.empty else set(key["key_interaction_id"])
    else:
        key_bits = set(key["key_interaction_id"]) if not key.empty else set()
    rows: list[dict[str, Any]] = []
    for fp in docked_fps.itertuples(index=False):
        pose = score_meta.loc[fp.pose_id]
        native_receptor = str(fp.ligand_id).rsplit("_", 1)[0].lower()
        native_ref = native_by_receptor.get(native_receptor)
        pose_bits = _bits(fp.fingerprint_sparse_json)
        native_bits = _bits(native_ref.fingerprint_sparse_json) if native_ref is not None else set()
        native_t, native_rec, native_prec, native_f1, missing, extra, recovered = _metrics(pose_bits, native_bits, key_bits)
        con_t, con_rec, con_prec, con_f1, con_missing, con_extra, con_recovered = _metrics(pose_bits, consensus_bits, key_bits)
        recovery_label = _label(native_rec if native_rec is not None else con_rec, native_t if native_t is not None else con_t, config)
        recovered_roles = "|".join(recovered | con_recovered)
        rows.append(
            {
                "pose_id": fp.pose_id,
                "docking_task_id": fp.docking_task_id,
                "ligand_id": fp.ligand_id,
                "target_receptor_id": fp.target_receptor_id,
                "task_type": fp.task_type,
                "docking_engine": fp.docking_engine,
                "pose_rank": fp.pose_rank,
                "rmsd_symmetry_corrected": pose["rmsd_symmetry_corrected"],
                "stage3_pose_label": pose["stage3_pose_label"],
                "sanity_status": pose["sanity_status"],
                "native_reference_available_flag": native_ref is not None,
                "native_reference_complex_id": None if native_ref is None else native_ref.complex_id,
                "binding_mode_reference_id": "native_consensus",
                "ifp_tanimoto_to_native": native_t,
                "ifp_tanimoto_to_consensus": con_t,
                "key_interaction_recall_native": native_rec,
                "key_interaction_precision_native": native_prec,
                "key_interaction_f1_native": native_f1,
                "key_interaction_recall_consensus": con_rec,
                "key_interaction_precision_consensus": con_prec,
                "key_interaction_f1_consensus": con_f1,
                "missing_key_interactions_json": json.dumps(sorted(missing or con_missing)),
                "extra_key_interactions_json": json.dumps(sorted(extra or con_extra)),
                "recovered_key_interactions_json": json.dumps(sorted(recovered or con_recovered)),
                "hinge_interaction_recovered_flag": "hinge:" in recovered_roles,
                "catalytic_lys_glu_region_consistent_flag": "catalytic_lys:" in recovered_roles or "catalytic_glu:" in recovered_roles,
                "gatekeeper_region_consistent_flag": "gatekeeper:" in recovered_roles,
                "dfg_region_consistent_flag": "dfg_region:" in recovered_roles,
                "interaction_recovery_label": recovery_label,
                "warnings_json": json.dumps([] if native_ref is not None else ["native_reference_unavailable"]),
            }
        )
    frame = pd.DataFrame(rows, columns=INTERACTION_RECOVERY_COLUMNS)
    write_table(paths["processed"] / "interaction_recovery.parquet", frame)
    write_table(paths["processed"] / "interaction_recovery.csv", frame)
    return frame
