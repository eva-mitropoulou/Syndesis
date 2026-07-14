"""Run controlled graph-preservation checks for strict pose reconstruction."""
from __future__ import annotations

import tempfile
import sys
from pathlib import Path

import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from syndesis.stage5.pose_reconstruction import reconstruct_pose_sdf  # noqa: E402

OUTPUT = ROOT / "results" / "robustness" / "pose_graph_mapping_validation.csv"
CASES = {
    "neutral_bond_orders": "CC(=O)Nc1ccncc1",
    "formal_charge": "C[N+](C)(C)CC(=O)[O-]",
    "tetrahedral_stereochemistry": "N[C@@H](C)C(=O)O",
    "alkene_stereochemistry": "C/C=C\\C(=O)N",
}


def write_pdbqt(path: Path, atoms: list[tuple[str, tuple[float, float, float]]]) -> None:
    lines = []
    for index, (element, (x, y, z)) in enumerate(atoms, start=1):
        lines.append(
            f"ATOM  {index:5d} {element + str(index):<4s} LIG A   1    "
            f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          {element:>2s}  {element}\n"
        )
    path.write_text("".join(lines) + "END\n")


def bond_signature(molecule: Chem.Mol) -> tuple[tuple[int, int, str], ...]:
    return tuple(sorted(
        (min(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()),
         max(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()), str(bond.GetBondType()))
        for bond in molecule.GetBonds()
    ))


def validate_case(case: str, smiles: str, root: Path) -> dict[str, object]:
    molecule = Chem.AddHs(Chem.MolFromSmiles(smiles))
    if AllChem.EmbedMolecule(molecule, randomSeed=807) != 0:
        raise RuntimeError(f"Embedding failed for validation case {case}")
    molecule = Chem.RemoveAllHs(molecule)
    template = root / f"{case}.sdf"
    writer = Chem.SDWriter(str(template))
    writer.write(molecule)
    writer.close()

    coordinates = molecule.GetConformer().GetPositions()
    elements = [atom.GetSymbol().upper() for atom in molecule.GetAtoms()]
    order = list(reversed(range(molecule.GetNumAtoms())))
    prepared_atoms = [(elements[index], tuple(coordinates[index])) for index in order]
    posed_atoms = [
        (element, (xyz[0] + 4.0, xyz[1] - 3.0, xyz[2] + 2.0))
        for element, xyz in prepared_atoms
    ]
    prepared = root / f"{case}_prepared.pdbqt"
    pose = root / f"{case}_pose.pdbqt"
    write_pdbqt(prepared, prepared_atoms)
    write_pdbqt(pose, posed_atoms)
    output = reconstruct_pose_sdf(pose, template, prepared, root / f"{case}_reconstructed.sdf")
    reconstructed = Chem.SDMolSupplier(str(output), removeHs=False)[0]
    reference_smiles = Chem.MolToSmiles(molecule, isomericSmiles=True)
    reconstructed_smiles = Chem.MolToSmiles(reconstructed, isomericSmiles=True)
    checks = {
        "atom_count_preserved": reconstructed.GetNumAtoms() == molecule.GetNumAtoms(),
        "graph_isomorphic": Chem.MolToSmiles(reconstructed, isomericSmiles=False)
        == Chem.MolToSmiles(molecule, isomericSmiles=False),
        "formal_charge_preserved": Chem.GetFormalCharge(reconstructed) == Chem.GetFormalCharge(molecule),
        "bond_orders_preserved": bond_signature(reconstructed) == bond_signature(molecule),
        "stereochemistry_preserved": reconstructed_smiles == reference_smiles,
        "coordinate_mapping_preserved": all(
            abs(reconstructed.GetConformer().GetAtomPosition(index).x - (coordinates[index][0] + 4.0)) <= 0.001
            and abs(reconstructed.GetConformer().GetAtomPosition(index).y - (coordinates[index][1] - 3.0)) <= 0.001
            and abs(reconstructed.GetConformer().GetAtomPosition(index).z - (coordinates[index][2] + 2.0)) <= 0.001
            for index in range(molecule.GetNumAtoms())
        ),
    }
    return {
        "case": case,
        "input_smiles": smiles,
        **{key: int(value) for key, value in checks.items()},
        "all_checks_pass": int(all(checks.values())),
    }
def main() -> int:
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory() as directory:
        rows = [validate_case(case, smiles, Path(directory)) for case, smiles in CASES.items()]
    report = pd.DataFrame(rows)
    if not report["all_checks_pass"].all():
        raise RuntimeError("At least one controlled graph-preservation check failed")
    report.to_csv(OUTPUT, index=False)
    print({"cases": len(report), "all_checks_pass": bool(report["all_checks_pass"].all())})
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
