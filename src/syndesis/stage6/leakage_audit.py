from __future__ import annotations

import pandas as pd

from syndesis.common.io import write_table
from syndesis.stage6 import schemas


def classify_feature(name: str, config: dict) -> tuple[str, bool, bool, str, str, str]:
    """Classify a feature column for training/deployment eligibility.

    This is a DEFAULT-DENY audit: a column is trainable only if it matches an
    explicitly-recognised deployment-safe family AND is not on the forbidden
    list. Unknown columns are quarantined (not trainable) rather than silently
    whitelisted, so that new label-derived or native-reference columns cannot
    leak into the model just because their name happens to contain a benign
    substring.
    """
    lowered = name.lower()
    metadata = set(schemas.IDENTIFIER_COLUMNS) | {"fingerprint_sparse_json", "fingerprint_bitstring", "warnings_json", "parse_warnings_json", "sanity_warnings_json"}
    forbidden_patterns = [str(x).lower() for x in config["features"].get("forbidden_patterns", [])]

    # 1. Hard leakage: explicit forbidden set, configured patterns, or any field
    #    that measures agreement with the native/consensus reference or the
    #    retrospective label (these define the training target).
    leakage_substrings = (
        "recovery", "recovered", "recall", "precision", "_f1", "tanimoto",
        "native", "consensus", "rmsd", "final_pose_label", "stage3_pose_label",
        "relevance", "binding_mode_reference", "cluster_label",
        "cluster_compatibility", "state_match",
    )
    if name in schemas.FORBIDDEN_FEATURES or any(pattern in lowered for pattern in forbidden_patterns):
        return "label_or_native_diagnostic", False, False, "high", "Native pose or retrospective label field.", "drop_for_training"
    if name.endswith("_json"):
        return "metadata", False, False, "medium", "Serialized list/blob retained outside the model matrix.", "metadata_only"
    if name in metadata or name.endswith("_id") or name.endswith("_file") or name.endswith("_path"):
        return "metadata", False, False, "medium", "Identifier or path retained outside the model matrix.", "metadata_only"
    if any(substr in lowered for substr in leakage_substrings):
        return "label_or_native_diagnostic", False, False, "high", "Encodes agreement with native/consensus reference or the label; target leakage.", "drop_for_training"

    # 2. Explicit deployment-safe allow-list (things computable for a brand-new
    #    docked pose with no knowledge of the native answer).
    if lowered.startswith("ifp_bit_"):
        return "prolif_interaction", True, True, "low", "Raw per-pose ProLIF interaction bit (describes the pose itself).", "train"
    if "cnn" in lowered or "vina" in lowered or "vinardo" in lowered or "gnina" in lowered or lowered.endswith("_score") or "affinity" in lowered or "score_disagreement" in lowered or "rank" in lowered:
        return "docking_rescoring", True, True, "low", "Deployment-safe docking or GNINA score feature.", "train"
    if "sanity" in lowered or "clash" in lowered or "strain" in lowered or lowered.endswith("_pocket") or "rg_" in lowered or "radius_of_gyration" in lowered:
        return "pose_sanity", True, True, "low", "Deployment-safe pose sanity/geometry feature.", "train"
    if "ligand_" in lowered or lowered.startswith("mol") or "tpsa" in lowered or "logp" in lowered or "hbd" in lowered or "hba" in lowered or "molecular_weight" in lowered or "rotatable" in lowered or "aromatic" in lowered:
        return "ligand_descriptor", True, True, "low", "Deployment-safe ligand descriptor.", "train"
    if "dfg_state" in lowered or "chelix" in lowered or "hrd" in lowered or "saltbridge" in lowered or "activation_loop" in lowered or "mutation" in lowered or "active_site_completeness" in lowered or "docking_box" in lowered or "receptor_state" in lowered or "suggested_docking_box" in lowered:
        return "receptor_state", True, True, "low", "Deployment-safe receptor-state metadata.", "train"

    # 3. Default deny: quarantine unknown columns instead of training on them.
    return "unrecognized_quarantined", False, False, "medium", "Unrecognized column quarantined by default-deny leakage policy.", "quarantine_review"


def audit_pose_features(features: pd.DataFrame, config: dict, paths: dict) -> pd.DataFrame:
    rows = []
    for col in features.columns:
        group, train, deploy, risk, reason, action = classify_feature(col, config)
        rows.append(
            {
                "feature_name": col,
                "feature_group": group,
                "allowed_for_training": bool(train),
                "allowed_for_deployment": bool(deploy),
                "leakage_risk": risk,
                "reason": reason,
                "action": action,
            }
        )
    audit = pd.DataFrame(rows, columns=schemas.LEAKAGE_AUDIT_COLUMNS)
    write_table(paths["processed"] / "feature_leakage_audit.parquet", audit)
    write_table(paths["processed"] / "feature_leakage_audit.csv", audit)
    return audit


def assert_no_leakage(training_columns: list[str]) -> None:
    bad = sorted(set(training_columns) & schemas.FORBIDDEN_FEATURES)
    if bad:
        raise RuntimeError(f"Forbidden leakage features reached training matrix: {bad}")
    banned_substrings = (
        "rmsd", "native_reference", "final_pose_label", "stage3_pose_label",
        "recovery", "recovered", "recall", "precision", "_f1", "tanimoto",
        "consensus", "relevance", "cluster_compatibility", "state_match",
        "binding_mode_reference",
    )
    for col in training_columns:
        lowered = col.lower()
        if any(substr in lowered for substr in banned_substrings):
            raise RuntimeError(f"Forbidden leakage feature reached training matrix: {col}")
