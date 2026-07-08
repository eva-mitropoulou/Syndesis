from __future__ import annotations

from egfr_dockingforge.stage7.schemas import ALLOWED_SOURCES


def test_source_tags_keep_vendor_and_known_separate() -> None:
    assert "zinc_vendor" in ALLOWED_SOURCES
    assert "chembl_known_ligand" in ALLOWED_SOURCES
    assert "generated_analog" in ALLOWED_SOURCES
