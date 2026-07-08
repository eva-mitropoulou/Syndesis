from __future__ import annotations

import pandas as pd

from egfr_dockingforge.stage5.prolif_engine import (
    compute_interactions,
    fingerprint_from_interactions,
    interaction_config_hash,
    prepare_ligand_for_prolif,
    prepare_protein_for_prolif,
)


def test_stage5_prolif_engine_returns_interactions() -> None:
    protein = prepare_protein_for_prolif(
        "data/processed/stage1/reference_complexes/1m17_a_aq4_999/receptor_clean.pdb",
        "data/processed/stage5/test_prolif_receptors",
    )
    ligand = prepare_ligand_for_prolif(
        "data/processed/stage1/reference_complexes/1m17_a_aq4_999/native_ligand.pdb",
        "data/processed/stage5/test_prolif_ligands",
    )
    residue_map = pd.read_csv("data/processed/stage2/pocket_residue_mapping.csv")
    residue_map = residue_map[residue_map["receptor_id"].eq("1m17_a_aq4_999")]
    config = {"interactions": {"enabled_interactions": ["Hydrophobic", "HBAcceptor", "HBDonor", "VdWContact"], "distance_cutoffs": {"hydrogen_bond": 3.6, "vdw": 4.0, "hydrophobic": 4.5}, "require_hydrogens": True}}
    interactions, meta = compute_interactions(protein, ligand, residue_map, config)
    assert not interactions.empty
    assert meta["interaction_engine"] == "prolif"
    assert meta["interaction_config_hash"] == interaction_config_hash(config)
    assert meta["warnings"] == []
    _bitstring, _sparse, bits = fingerprint_from_interactions(interactions)
    assert bits


def test_stage5_ligand_preparation_uses_unique_source_paths() -> None:
    first = prepare_ligand_for_prolif(
        "data/processed/stage1/reference_complexes/1m17_a_aq4_999/native_ligand.pdb",
        "data/processed/stage5/test_prolif_ligands",
    )
    second = prepare_ligand_for_prolif(
        "data/processed/stage1/reference_complexes/4hjo_a_aq4_1001/native_ligand.pdb",
        "data/processed/stage5/test_prolif_ligands",
    )
    assert first.name != second.name
    assert first.exists()
    assert second.exists()
