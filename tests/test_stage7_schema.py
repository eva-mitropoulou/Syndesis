from __future__ import annotations

from egfr_dockingforge.stage7 import schemas


def test_stage7_schema_columns() -> None:
    assert "molecule_id" in schemas.MASTER_COLUMNS
    assert "prepared_ligand_id" in schemas.STAGE8_COLUMNS
    assert "zinc_vendor" in schemas.ALLOWED_SOURCES
