from __future__ import annotations

from egfr_dockingforge.stage2.pocket_mapping import mapping_frame


def test_missing_pocket_mapping_records_warning() -> None:
    rows = [
        {
            "receptor_id": "r1",
            "pdb_id": "TEST",
            "auth_asym_id": "A",
            "uniprot_residue_number": 745,
            "auth_seq_id": None,
            "residue_name": None,
            "klifs_position": None,
            "atom_names_present": "",
            "residue_complete_flag": False,
            "used_for_alignment_flag": False,
            "warning": "Residue missing under auth_seq_id fallback mapping.",
        }
    ]
    frame = mapping_frame(rows)
    assert frame.loc[0, "warning"]
    assert frame.loc[0, "residue_complete_flag"] is False or frame.loc[0, "residue_complete_flag"] == False

