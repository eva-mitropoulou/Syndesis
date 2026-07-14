from pathlib import Path

import pytest
from rdkit import Chem
from rdkit.Chem import AllChem

from syndesis.stage5.pose_reconstruction import reconstruct_pose_sdf


def _write_pdbqt(path: Path, atoms: list[tuple[str, tuple[float, float, float]]]) -> None:
    lines = []
    for index, (element, (x, y, z)) in enumerate(atoms, start=1):
        name = f"{element}{index}"
        lines.append(
            f"ATOM  {index:5d} {name:<4s} LIG A   1    "
            f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          {element:>2s}  {element}\n"
        )
    path.write_text("".join(lines) + "END\n")


def test_reconstruct_pose_preserves_template_graph_and_maps_reordered_atoms(tmp_path: Path) -> None:
    molecule = Chem.AddHs(Chem.MolFromSmiles("CCO"))
    AllChem.EmbedMolecule(molecule, randomSeed=7)
    molecule = Chem.RemoveHs(molecule)
    template = tmp_path / "lig.sdf"
    writer = Chem.SDWriter(str(template))
    writer.write(molecule)
    writer.close()

    coordinates = molecule.GetConformer().GetPositions()
    elements = [atom.GetSymbol().upper() for atom in molecule.GetAtoms()]
    order = [2, 0, 1]
    prepared_atoms = [(elements[index], tuple(coordinates[index])) for index in order]
    pose_atoms = [(element, (xyz[0] + 10.0, xyz[1] - 2.0, xyz[2] + 1.0)) for element, xyz in prepared_atoms]
    prepared = tmp_path / "lig.pdbqt"
    pose = tmp_path / "pose.pdbqt"
    _write_pdbqt(prepared, prepared_atoms)
    _write_pdbqt(pose, pose_atoms)

    output = reconstruct_pose_sdf(pose, template, prepared, tmp_path / "pose.sdf")
    reconstructed = Chem.SDMolSupplier(str(output), removeHs=True)[0]
    assert reconstructed is not None
    assert Chem.MolToSmiles(reconstructed) == Chem.MolToSmiles(molecule)
    result = reconstructed.GetConformer().GetPositions()
    for index in range(molecule.GetNumAtoms()):
        assert result[index] == pytest.approx(coordinates[index] + [10.0, -2.0, 1.0], abs=0.002)


def test_reconstruct_pose_rejects_unmappable_prepared_coordinates(tmp_path: Path) -> None:
    molecule = Chem.MolFromSmiles("CO")
    AllChem.Compute2DCoords(molecule)
    template = tmp_path / "lig.sdf"
    writer = Chem.SDWriter(str(template))
    writer.write(molecule)
    writer.close()
    prepared = tmp_path / "lig.pdbqt"
    pose = tmp_path / "pose.pdbqt"
    atoms = [("C", (100.0, 0.0, 0.0)), ("O", (101.0, 0.0, 0.0))]
    _write_pdbqt(prepared, atoms)
    _write_pdbqt(pose, atoms)
    with pytest.raises(RuntimeError, match="mapping exceeds"):
        reconstruct_pose_sdf(pose, template, prepared, tmp_path / "pose.sdf")


def test_reconstruct_pose_makes_stereo_defining_hydrogen_implicit(tmp_path: Path) -> None:
    molecule = Chem.AddHs(Chem.MolFromSmiles("[H]/N=C(/C)\\C"))
    AllChem.EmbedMolecule(molecule, randomSeed=807)
    template = tmp_path / "template.sdf"
    writer = Chem.SDWriter(str(template))
    writer.write(molecule)
    writer.close()

    heavy = Chem.RemoveAllHs(molecule)
    coordinates = heavy.GetConformer().GetPositions()
    atoms = [
        (atom.GetSymbol().upper(), tuple(coordinates[index]))
        for index, atom in enumerate(heavy.GetAtoms())
    ]
    prepared = tmp_path / "prepared.pdbqt"
    pose = tmp_path / "pose.pdbqt"
    _write_pdbqt(prepared, atoms)
    _write_pdbqt(
        pose,
        [(element, (xyz[0] + 2.0, xyz[1] + 2.0, xyz[2] + 2.0)) for element, xyz in atoms],
    )

    output = reconstruct_pose_sdf(pose, template, prepared, tmp_path / "pose.sdf")
    reconstructed = Chem.SDMolSupplier(str(output), removeHs=False)[0]

    assert reconstructed.GetNumAtoms() == heavy.GetNumAtoms()
    assert Chem.MolToSmiles(reconstructed) == Chem.MolToSmiles(heavy)
    assert sum(atom.GetTotalNumHs() for atom in reconstructed.GetAtoms()) > 0


@pytest.mark.parametrize(
    "smiles",
    [
        "C[N+](C)(C)CC(=O)[O-]",
        "N[C@@H](C)C(=O)O",
        "C/C=C\\C(=O)N",
    ],
)
def test_reconstruct_pose_preserves_charge_bonds_and_stereochemistry(
    tmp_path: Path, smiles: str
) -> None:
    molecule = Chem.AddHs(Chem.MolFromSmiles(smiles))
    assert AllChem.EmbedMolecule(molecule, randomSeed=807) == 0
    molecule = Chem.RemoveAllHs(molecule)
    template = tmp_path / "template.sdf"
    writer = Chem.SDWriter(str(template))
    writer.write(molecule)
    writer.close()
    coordinates = molecule.GetConformer().GetPositions()
    elements = [atom.GetSymbol().upper() for atom in molecule.GetAtoms()]
    order = list(reversed(range(molecule.GetNumAtoms())))
    prepared_atoms = [(elements[index], tuple(coordinates[index])) for index in order]
    prepared = tmp_path / "prepared.pdbqt"
    pose = tmp_path / "pose.pdbqt"
    _write_pdbqt(prepared, prepared_atoms)
    _write_pdbqt(
        pose,
        [(element, (xyz[0] + 1.0, xyz[1] + 2.0, xyz[2] + 3.0)) for element, xyz in prepared_atoms],
    )

    output = reconstruct_pose_sdf(pose, template, prepared, tmp_path / "pose.sdf")
    reconstructed = Chem.SDMolSupplier(str(output), removeHs=False)[0]

    assert Chem.MolToSmiles(reconstructed, isomericSmiles=True) == Chem.MolToSmiles(
        molecule, isomericSmiles=True
    )
    assert Chem.GetFormalCharge(reconstructed) == Chem.GetFormalCharge(molecule)
    assert [str(bond.GetBondType()) for bond in reconstructed.GetBonds()] == [
        str(bond.GetBondType()) for bond in molecule.GetBonds()
    ]
