from __future__ import annotations

from Bio.PDB import PDBParser

from egfr_dockingforge.stage2.pocket_mapping import resolve_uniprot_residue, residue_by_auth_seq


def test_resolve_uniprot_residue_handles_egfr_auth_seq_offset() -> None:
    structure = PDBParser(QUIET=True).get_structure(
        "1m17",
        "data/processed/stage1/reference_complexes/1m17_a_aq4_999/receptor_clean.pdb",
    )
    chain = next(structure.get_chains())
    residues = residue_by_auth_seq(chain)

    residue, auth_seq_id, mapping_method = resolve_uniprot_residue(residues, 745)

    assert residue is not None
    assert residue.get_resname().strip() == "LYS"
    assert auth_seq_id == 721
    assert mapping_method == "egfr_sequence_offset_-24"
