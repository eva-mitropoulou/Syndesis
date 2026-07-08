import pandas as pd

from egfr_dockingforge.stage9.analog_validation import validate_analog_batch


def test_invalid_smiles_rejected_and_scope_flags_stored(tmp_path):
    candidates = pd.DataFrame(
        [
            {"analog_id": "a1", "standard_smiles": "not_smiles"},
            {"analog_id": "a2", "standard_smiles": "C=CC(=O)Nc1ccccc1"},
            {"analog_id": "a3", "standard_smiles": "CCOc1ccc(Nc2ncnc3ccccc23)cc1"},
        ]
    )
    out = validate_analog_batch(candidates, {}, {"processed": tmp_path})
    assert out.loc[out["analog_id"].eq("a1"), "rejection_reason"].item() == "invalid_smiles"
    assert bool(out.loc[out["analog_id"].eq("a2"), "covalent_warhead_flag"].item())
    assert bool(out.loc[out["analog_id"].eq("a3"), "hard_scope_pass"].item())
