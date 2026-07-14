from __future__ import annotations

import numpy as np
import pytest

from syndesis.stage3 import rmsd
from syndesis.stage3.rmsd import mapped_rmsd


def test_mapped_rmsd_refuses_mismatched_atom_counts() -> None:
    with pytest.raises(ValueError):
        mapped_rmsd(np.zeros((2, 3)), np.zeros((3, 3)))


def test_mapped_rmsd_zero_for_identical_coordinates() -> None:
    coords = np.array([[0, 0, 0], [1, 0, 0]], dtype=float)
    assert mapped_rmsd(coords, coords) == 0.0


def test_symmetry_rmsd_rejects_missing_template_graph(monkeypatch) -> None:
    monkeypatch.setattr(rmsd, "_load_mol_with_coords", lambda *_args: None)
    with pytest.raises(RuntimeError, match="template-mapped molecular graphs"):
        rmsd.symmetry_corrected_rmsd("pose.pdbqt", "reference.pdb", None)
