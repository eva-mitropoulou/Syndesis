from __future__ import annotations

from syndesis.stage1.complex_classification import smiles_has_warhead
from syndesis.stage1.rcsb_client import search_pdb_ids
from syndesis.stage0.scope_schema import load_yaml_mapping


def test_covalent_warhead_smiles_are_detected() -> None:
    flag, reason = smiles_has_warhead("C=CC(=O)Nc1ccc2ncnc(N)c2c1")
    assert flag is True
    assert reason


def test_stage1_config_excludes_covalent_and_cofactor_ligands() -> None:
    config = load_yaml_mapping("configs/stage1_cocrystal_benchmark.yaml")
    known_covalent = set(config["filters"]["known_covalent_ligand_comp_ids"])
    excluded = set(config["filters"]["excluded_ligand_comp_ids"])
    assert "4ZAU" in config["controls"]["expected_scope_exclusion_controls"]
    assert known_covalent
    assert {"ATP", "ADP", "ANP"}.issubset(excluded)


def test_reference_controls_are_in_candidate_ids() -> None:
    config = load_yaml_mapping("configs/stage1_cocrystal_benchmark.yaml")
    ids = set(search_pdb_ids(config))
    assert {"1M17", "1XKK", "4ZAU"}.issubset(ids)

