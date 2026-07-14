from __future__ import annotations

from dataclasses import dataclass
from typing import Any

from Bio.PDB import NeighborSearch
from Bio.PDB.Polypeptide import is_aa

from syndesis.stage1.mmcif_parser import WATER_RESNAMES


@dataclass(frozen=True)
class LigandInstance:
    comp_id: str
    instance_id: str
    auth_asym_id: str
    label_asym_id: str | None
    residue_id: tuple[str, int, str]
    heavy_atom_count: int
    occupancy_min: float | None
    occupancy_mean: float | None
    bfactor_mean: float | None
    altloc_flag: bool
    min_distance_to_chain: float | None
    pocket_bfactor_mean: float | None


def atom_element(atom: Any) -> str:
    return (atom.element or atom.get_name()[0]).strip().upper()


def heavy_atoms(residue: Any) -> list[Any]:
    return [atom for atom in residue.get_atoms() if atom_element(atom) != "H"]


def ligand_stats(residue: Any) -> tuple[int, float | None, float | None, float | None, bool]:
    atoms = heavy_atoms(residue)
    occupancies = [atom.get_occupancy() for atom in atoms if atom.get_occupancy() is not None]
    bfactors = [atom.get_bfactor() for atom in atoms if atom.get_bfactor() is not None]
    altloc = any(str(atom.get_altloc()).strip() not in {"", " "} for atom in atoms)
    occ_min = min(occupancies) if occupancies else None
    occ_mean = sum(occupancies) / len(occupancies) if occupancies else None
    b_mean = sum(bfactors) / len(bfactors) if bfactors else None
    return len(atoms), occ_min, occ_mean, b_mean, altloc


def protein_atoms(chain: Any) -> list[Any]:
    return [
        atom
        for residue in chain.get_residues()
        if is_aa(residue, standard=True)
        for atom in residue.get_atoms()
        if atom_element(atom) != "H"
    ]


def distance_and_pocket_bfactor(residue: Any, chain: Any, cutoff: float = 5.0) -> tuple[float | None, float | None]:
    atoms = heavy_atoms(residue)
    chain_atoms = protein_atoms(chain)
    if not atoms or not chain_atoms:
        return None, None
    search = NeighborSearch(chain_atoms)
    distances: list[float] = []
    pocket_b: list[float] = []
    for atom in atoms:
        neighbors = search.search(atom.coord, cutoff, level="A")
        for neighbor in neighbors:
            distances.append(atom - neighbor)
            pocket_b.append(neighbor.get_bfactor())
    min_distance = min(distances) if distances else None
    pocket_mean = sum(pocket_b) / len(pocket_b) if pocket_b else None
    return min_distance, pocket_mean


def residue_instance_id(residue: Any) -> str:
    het, seq, icode = residue.id
    suffix = str(icode).strip()
    return f"{seq}{suffix}" if suffix else str(seq)


def is_small_molecule_residue(residue: Any, excluded_comp_ids: set[str]) -> bool:
    comp_id = residue.get_resname().strip().upper()
    het_flag = str(residue.id[0]).strip()
    if not het_flag:
        return False
    if comp_id in WATER_RESNAMES or comp_id in excluded_comp_ids:
        return False
    if is_aa(residue, standard=True):
        return False
    atoms = heavy_atoms(residue)
    if not atoms:
        return False
    return any(atom_element(atom) == "C" for atom in atoms)


def find_ligands_for_chain(chain: Any, config: dict[str, Any]) -> list[LigandInstance]:
    excluded = {str(value).upper() for value in config["filters"].get("excluded_ligand_comp_ids", [])}
    cutoff = float(config["filters"].get("max_ligand_to_pocket_distance_angstrom", 5.0))
    ligands: list[LigandInstance] = []
    for residue in chain.get_residues():
        if not is_small_molecule_residue(residue, excluded):
            continue
        heavy_count, occ_min, occ_mean, b_mean, altloc = ligand_stats(residue)
        min_distance, pocket_b = distance_and_pocket_bfactor(residue, chain, cutoff=cutoff)
        ligands.append(
            LigandInstance(
                comp_id=residue.get_resname().strip().upper(),
                instance_id=residue_instance_id(residue),
                auth_asym_id=chain.id,
                label_asym_id=None,
                residue_id=residue.id,
                heavy_atom_count=heavy_count,
                occupancy_min=occ_min,
                occupancy_mean=occ_mean,
                bfactor_mean=b_mean,
                altloc_flag=altloc,
                min_distance_to_chain=min_distance,
                pocket_bfactor_mean=pocket_b,
            )
        )
    ligands.sort(key=lambda item: (item.min_distance_to_chain is not None, item.heavy_atom_count), reverse=True)
    return ligands

