from __future__ import annotations

from pathlib import Path
from typing import Any

from Bio.PDB import PDBIO, Select
from Bio.PDB.Polypeptide import is_aa

from syndesis.common.io import ensure_dir


class ChainSelect(Select):
    def __init__(self, chain_id: str) -> None:
        self.chain_id = chain_id

    def accept_chain(self, chain: Any) -> bool:
        return chain.id == self.chain_id


class ReceptorSelect(ChainSelect):
    def accept_residue(self, residue: Any) -> bool:
        return bool(is_aa(residue, standard=True))


class LigandSelect(ChainSelect):
    def __init__(self, chain_id: str, residue_id: tuple[str, int, str]) -> None:
        super().__init__(chain_id)
        self.residue_id = residue_id

    def accept_residue(self, residue: Any) -> bool:
        return residue.id == self.residue_id


class ComplexSelect(ChainSelect):
    def __init__(self, chain_id: str, residue_id: tuple[str, int, str]) -> None:
        super().__init__(chain_id)
        self.residue_id = residue_id

    def accept_residue(self, residue: Any) -> bool:
        return bool(is_aa(residue, standard=True)) or residue.id == self.residue_id


class PocketWaterSelect(ChainSelect):
    def accept_residue(self, residue: Any) -> bool:
        return residue.get_resname().strip().upper() in {"HOH", "WAT", "DOD", "H2O"}


def save_selection(structure: Any, path: Path, selector: Select) -> Path:
    ensure_dir(path.parent)
    io = PDBIO()
    io.set_structure(structure)
    io.save(str(path), selector)
    return path


def export_reference_files(
    structure: Any,
    complex_id: str,
    chain_id: str,
    ligand_residue_id: tuple[str, int, str],
    out_root: Path,
) -> dict[str, str]:
    out_dir = ensure_dir(out_root / complex_id.lower())
    native_complex = save_selection(structure, out_dir / "native_complex.pdb", ComplexSelect(chain_id, ligand_residue_id))
    receptor_clean = save_selection(structure, out_dir / "receptor_clean.pdb", ReceptorSelect(chain_id))
    native_ligand = save_selection(structure, out_dir / "native_ligand.pdb", LigandSelect(chain_id, ligand_residue_id))
    pocket_waters = save_selection(structure, out_dir / "pocket_waters.pdb", PocketWaterSelect(chain_id))
    return {
        "native_complex_path": str(native_complex),
        "native_ligand_sdf_path": str(native_ligand),
        "receptor_clean_path": str(receptor_clean),
        "pocket_waters_path": str(pocket_waters),
    }

