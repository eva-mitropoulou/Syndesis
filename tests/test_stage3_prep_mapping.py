from __future__ import annotations

from egfr_dockingforge.stage3.task_matrix import ca_atoms
from egfr_dockingforge.stage3.rmsd import mapped_rmsd, heavy_atom_coords


def test_stage3_ca_atoms_use_corrected_egfr_residue_mapping() -> None:
    atoms = ca_atoms(
        "data/processed/stage1/reference_complexes/1m17_a_aq4_999/receptor_clean.pdb",
        [745, 762, 790, 793, 797, 855, 856, 857],
    )

    assert set(atoms) == {745, 762, 790, 793, 797, 855, 856, 857}


def test_stage3_rmsd_uses_pdbqt_reference_atom_order() -> None:
    reference = "data/processed/stage3/docking_ligands/5cav_a_4zq_1101_4zq.pdbqt"
    pose = "data/processed/stage3/docking_outputs/poses/run__5cav_a_4zq_1101_4zq__5cav_a_4zq_1101__unidock__seed13__pose01.pdbqt"

    assert mapped_rmsd(heavy_atom_coords(reference), heavy_atom_coords(pose)) < 1.0
