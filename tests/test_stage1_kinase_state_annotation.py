from __future__ import annotations

from egfr_dockingforge.stage1.kincore_metadata import infer_kincore_metadata


def test_stage1_uses_klifs_dfg_ac_helix_for_active_state() -> None:
    metadata, warnings = infer_kincore_metadata(
        "5CAV",
        "A",
        {},
        {"klifs_dfg_state": "in", "klifs_ac_helix_state": "in"},
    )
    assert warnings == []
    assert metadata["kincore_state"] == "active-like"
    assert metadata["dfg_state"] == "DFGin"
    assert metadata["chelix_state"] == "C-helix-in"


def test_stage1_uses_klifs_dfg_ac_helix_for_inactive_state() -> None:
    metadata, warnings = infer_kincore_metadata(
        "6DUK",
        "A",
        {},
        {"klifs_dfg_state": "in", "klifs_ac_helix_state": "out"},
    )
    assert warnings == []
    assert metadata["kincore_state"] == "inactive-like"
    assert metadata["dfg_state"] == "DFGin"
    assert metadata["chelix_state"] == "C-helix-out"


def test_stage1_uses_klifs_dfgout_for_typeii_like_state() -> None:
    metadata, warnings = infer_kincore_metadata(
        "5HG5",
        "A",
        {},
        {"klifs_dfg_state": "out", "klifs_ac_helix_state": "out"},
    )
    assert warnings == []
    assert metadata["kincore_state"] == "inactive-like"
    assert metadata["dfg_state"] == "DFGout"
