"""Strictly transfer docked coordinates onto a prepared ligand molecular graph."""
from __future__ import annotations

from pathlib import Path

from rdkit import Chem
from rdkit.Geometry import Point3D

from egfr_dockingforge.stage5.prolif_engine import _autodock_element


def _heavy_atoms(path: str | Path) -> list[tuple[str, tuple[float, float, float]]]:
    atoms = []
    for raw in Path(path).read_text(encoding="utf-8", errors="ignore").splitlines():
        if not raw.startswith(("ATOM", "HETATM")):
            continue
        element = _autodock_element(raw)
        if element == "H":
            continue
        atoms.append((element, (float(raw[30:38]), float(raw[38:46]), float(raw[46:54]))))
    return atoms


def _heavy_atom_template(path: str | Path) -> Chem.Mol:
    """Load the prepared graph with every explicit hydrogen made implicit.

    RDKit retains hydrogens that define double-bond stereochemistry when an SDF
    is read with ``removeHs=True``. OpenBabel writes those hydrogens implicitly
    in PDBQT, so ``RemoveAllHs`` is required for an exact heavy-atom comparison.
    Formal charges, tautomerism, and hydrogen counts remain encoded on the heavy
    atoms.
    """
    supplier = Chem.SDMolSupplier(str(path), removeHs=False)
    molecule = supplier[0] if supplier and len(supplier) else None
    if molecule is None:
        raise RuntimeError(f"Unreadable prepared SDF template: {path}")
    return Chem.RemoveAllHs(molecule)


def reconstruct_pose_sdf(
    pose_file: str | Path,
    template_sdf: str | Path,
    prepared_pdbqt: str | Path,
    output_sdf: str | Path,
    *,
    mapping_tolerance_angstrom: float = 0.05,
) -> Path:
    """Return a valid SDF with the template graph and unchanged docked coordinates.

    Open Babel can reorder atoms while converting SDF to PDBQT. The template-to-
    PDBQT map is therefore recovered from the pre-docking conformer coordinates,
    with exact element matching and a strict distance tolerance. Docking preserves
    PDBQT atom order, allowing the corresponding pose coordinates to be transferred.
    """
    template = _heavy_atom_template(template_sdf)
    prepared = _heavy_atoms(prepared_pdbqt)
    posed = _heavy_atoms(pose_file)
    if len(prepared) != len(posed) or len(prepared) != template.GetNumAtoms():
        raise RuntimeError(
            f"Atom-count mismatch: template={template.GetNumAtoms()}, "
            f"prepared={len(prepared)}, posed={len(posed)}"
        )
    if [atom[0] for atom in prepared] != [atom[0] for atom in posed]:
        raise RuntimeError("Prepared/posed atom order or element sequence changed")

    template_coords = template.GetConformer().GetPositions()
    mapping: dict[int, int] = {}
    unused = set(range(len(prepared)))
    for template_index, atom in enumerate(template.GetAtoms()):
        element = atom.GetSymbol().upper()
        candidates = [index for index in unused if prepared[index][0] == element]
        if not candidates:
            raise RuntimeError(f"No prepared {element} atom for template atom {template_index}")
        distances = [
            float(sum((template_coords[template_index][axis] - prepared[index][1][axis]) ** 2 for axis in range(3)) ** 0.5)
            for index in candidates
        ]
        offset = min(range(len(candidates)), key=distances.__getitem__)
        best = candidates[offset]
        if distances[offset] > mapping_tolerance_angstrom:
            raise RuntimeError(
                f"Template/prepared mapping exceeds {mapping_tolerance_angstrom:.3f} A: "
                f"atom {template_index}, distance {distances[offset]:.4f} A"
            )
        mapping[template_index] = best
        unused.remove(best)
    if unused:
        raise RuntimeError(f"Incomplete atom mapping: {sorted(unused)}")

    molecule = Chem.Mol(template)
    conformer = Chem.Conformer(molecule.GetNumAtoms())
    for template_index, prepared_index in mapping.items():
        conformer.SetAtomPosition(template_index, Point3D(*posed[prepared_index][1]))
    molecule.RemoveAllConformers()
    molecule.AddConformer(conformer)
    target = Path(output_sdf)
    target.parent.mkdir(parents=True, exist_ok=True)
    writer = Chem.SDWriter(str(target))
    writer.write(molecule)
    writer.close()
    if not target.exists() or target.stat().st_size == 0:
        raise RuntimeError(f"Failed to write reconstructed pose SDF: {target}")
    return target
