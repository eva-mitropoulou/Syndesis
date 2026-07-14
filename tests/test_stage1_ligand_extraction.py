from __future__ import annotations

from types import SimpleNamespace

from syndesis.stage1.ligand_extraction import ligand_stats


class FakeAtom:
    element = "C"

    def __init__(self, name: str = "C1") -> None:
        self.name = name

    def get_occupancy(self) -> float:
        return 1.0

    def get_bfactor(self) -> float:
        return 20.0

    def get_altloc(self) -> str:
        return " "

    def get_name(self) -> str:
        return self.name


def test_ligand_atom_count_is_nonzero_for_extracted_residue() -> None:
    residue = SimpleNamespace(get_atoms=lambda: [FakeAtom(), FakeAtom("N1")])
    heavy_count, occ_min, occ_mean, bfactor_mean, altloc = ligand_stats(residue)
    assert heavy_count == 2
    assert occ_min == 1.0
    assert occ_mean == 1.0
    assert bfactor_mean == 20.0
    assert altloc is False

