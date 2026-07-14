from __future__ import annotations

import pandas as pd

from syndesis.stage5.native_atlas import build_key_interaction_map


def test_stage5_key_interactions_require_native_support() -> None:
    native = pd.DataFrame(
        [
            {"complex_id": "a", "residue_role": "hinge", "uniprot_residue_number": 793, "auth_seq_id": 793, "interaction_type": "HBAcceptor", "residue_name": "MET", "klifs_position": None},
            {"complex_id": "b", "residue_role": "hinge", "uniprot_residue_number": 793, "auth_seq_id": 793, "interaction_type": "HBAcceptor", "residue_name": "MET", "klifs_position": None},
            {"complex_id": "a", "residue_role": "pocket_residue", "uniprot_residue_number": 800, "auth_seq_id": 800, "interaction_type": "Hydrophobic", "residue_name": "LEU", "klifs_position": None},
        ]
    )
    config = {"key_interactions": {"key_interaction_threshold_native_frequency": 0.25, "minimum_native_count": 2, "key_interaction_manual_overrides": []}}
    key = build_key_interaction_map(native, config)
    assert "hinge:793:HBAcceptor" in set(key["key_interaction_id"])
    assert "pocket_residue:800:Hydrophobic" not in set(key["key_interaction_id"])

