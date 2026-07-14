from __future__ import annotations

import json
import pandas as pd


def select_subsets(master: pd.DataFrame, config: dict) -> pd.DataFrame:
    rows = []
    candidates = master[master["include_in_screening_library"]].copy()
    controls = master[master["screening_role"].isin(["native_pose_reference", "known_activity_reference"])]
    for _, row in controls.iterrows():
        rows.append({**_base(row), "screening_subset": "known_controls", "selection_reason": "known/native reference", "include_in_subset": True})
        rows.append({**_base(row), "screening_subset": "retrospective_known_ligands", "selection_reason": "known EGFR ligand", "include_in_subset": True})
    for subset, size in [("smoke_test_1k", config["diversity"]["smoke_test_size"]), ("dev_screen_10k_or_50k", config["diversity"]["dev_screen_size"]), ("main_screen_100k_or_250k", config["diversity"]["main_screen_size"])]:
        for _, row in candidates.head(size).iterrows():
            rows.append({**_base(row), "screening_subset": subset, "selection_reason": "deterministic source-aware selection", "include_in_subset": True})
    for bucket, subset in [("close_analog", "vendor_close_analogs"), ("scaffold_novel", "vendor_hard_novel")]:
        for _, row in candidates[candidates["novelty_bucket"].eq(bucket)].iterrows():
            rows.append({**_base(row), "screening_subset": subset, "selection_reason": bucket, "include_in_subset": True})
    for src, subset in [("zinc_vendor", "vendor_diverse"), ("generated_analog", "generated_analogs"), ("manual_analog", "manual_analogs")]:
        for _, row in candidates[candidates["source"].eq(src)].iterrows():
            rows.append({**_base(row), "screening_subset": subset, "selection_reason": src, "include_in_subset": True})
    return pd.DataFrame(rows)


def _base(row) -> dict:
    return {
        "molecule_id": row["molecule_id"], "source": row["source"], "subsource": row["subsource"],
        "novelty_bucket": row["novelty_bucket"], "scaffold_id": row["scaffold_id"], "cluster_id": row["scaffold_id"],
        "cluster_medoid_flag": True, "property_bin": "stage7_default", "diversity_score": 1.0,
        "warnings_json": json.dumps([]),
    }
