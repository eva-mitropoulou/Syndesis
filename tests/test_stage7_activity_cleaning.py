from __future__ import annotations

from egfr_dockingforge.stage7.activity_cleaning import p_activity_from_nm, to_nm


def test_units_convert_to_nm_and_pactivity() -> None:
    assert to_nm(1, "uM") == 1000
    assert round(p_activity_from_nm(1000), 3) == 6.0
