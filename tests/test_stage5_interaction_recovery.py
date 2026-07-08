from __future__ import annotations

from pathlib import Path

import pandas as pd

from egfr_dockingforge.stage5.interaction_recovery import compute_interaction_recovery


def test_stage5_interaction_recovery_calculates_metrics(tmp_path: Path) -> None:
    docked = pd.DataFrame([{"pose_id": "p1", "docking_task_id": "t1", "ligand_id": "native_lig", "target_receptor_id": "rec", "task_type": "redocking", "docking_engine": "unidock", "pose_rank": 1, "fingerprint_sparse_json": '["hinge:793:HBAcceptor"]'}])
    native = pd.DataFrame([{"complex_id": "native", "receptor_id": "native", "fingerprint_sparse_json": '["hinge:793:HBAcceptor","gatekeeper:790:Hydrophobic"]'}])
    key = pd.DataFrame([{"key_interaction_id": "hinge:793:HBAcceptor"}, {"key_interaction_id": "gatekeeper:790:Hydrophobic"}])
    inputs = {"pose_scores": pd.DataFrame([{"pose_id": "p1", "rmsd_symmetry_corrected": 1.2, "stage3_pose_label": "strict_native_like", "sanity_status": "pass"}])}
    config = {"recovery": {"high_key_recall_threshold": 0.75, "moderate_key_recall_threshold": 0.5, "high_tanimoto_threshold": 0.5}}
    rec = compute_interaction_recovery(docked, native, key, inputs, config, {"processed": tmp_path})
    assert rec.loc[0, "ifp_tanimoto_to_native"] == 0.5
    assert rec.loc[0, "key_interaction_recall_native"] == 0.5
    assert rec.loc[0, "interaction_recovery_label"] == "moderate_recovery"

