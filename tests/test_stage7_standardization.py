from __future__ import annotations

from syndesis.stage7.standardize import standardize_smiles


def test_standardization_invalid_smiles_fails_explicitly() -> None:
    row = standardize_smiles("not_a_smiles", "x", "manual_analog")
    assert row["standardization_status"] == "failed"
    assert "invalid_smiles" in row["warnings_json"]


def test_standardization_salt_keeps_largest_fragment() -> None:
    row = standardize_smiles("CCO.Cl", "x", "manual_analog")
    assert row["standardization_status"] == "success"
    assert row["salt_removed_flag"]
    assert row["inchi_key"]
