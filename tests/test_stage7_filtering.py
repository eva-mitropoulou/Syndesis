from __future__ import annotations

from egfr_dockingforge.stage7.filtering import flags


def test_covalent_warhead_excluded() -> None:
    config = {"filters": {"mw_min": 10, "mw_max": 1000, "heavy_atom_min": 1, "heavy_atom_max": 100, "clogp_min": -10, "clogp_max": 10, "tpsa_min": 0, "tpsa_max": 500, "hbd_max": 20, "hba_max": 30, "rotatable_bonds_max": 50, "formal_charge_min": -5, "formal_charge_max": 5, "exclude_macrocycles": True}}
    out = flags("C=CC(=O)N", "zinc_vendor", config)
    assert out["covalent_warhead_flag"]
    assert not out["include_in_screening_library"]
