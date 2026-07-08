from __future__ import annotations

import pandas as pd

from egfr_dockingforge.stage5.residue_mapping import normalize_interaction_key, residue_role


def test_stage5_residue_roles_include_key_egfr_residues() -> None:
    assert residue_role(745) == "catalytic_lys"
    assert residue_role(793) == "hinge"
    assert residue_role(857) == "dfg_region"


def test_stage5_interaction_key_uses_role_residue_and_type() -> None:
    key = normalize_interaction_key(pd.Series({"uniprot_residue_number": 793, "residue_role": "hinge", "interaction_type": "HBAcceptor"}))
    assert key == "hinge:793:HBAcceptor"

