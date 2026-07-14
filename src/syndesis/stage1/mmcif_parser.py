from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any

from Bio.PDB import MMCIFParser
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from Bio.PDB.Polypeptide import is_aa


WATER_RESNAMES = {"HOH", "WAT", "DOD", "H2O"}


@dataclass(frozen=True)
class ParsedMmcif:
    path: Path
    mmcif: dict[str, Any]
    structure: Any


def parse_mmcif(path: str | Path) -> ParsedMmcif:
    path = Path(path)
    parser = MMCIFParser(QUIET=True)
    return ParsedMmcif(path=path, mmcif=MMCIF2Dict(str(path)), structure=parser.get_structure(path.stem.upper(), str(path)))


def first_value(mmcif: dict[str, Any], key: str) -> Any | None:
    value = mmcif.get(key)
    if value is None:
        return None
    if isinstance(value, list):
        return value[0] if value else None
    return value


def as_list(value: Any) -> list[Any]:
    if value is None:
        return []
    return value if isinstance(value, list) else [value]


def float_or_none(value: Any) -> float | None:
    try:
        if value in {None, ".", "?"}:
            return None
        return float(value)
    except (TypeError, ValueError):
        return None


def entry_resolution(mmcif: dict[str, Any]) -> float | None:
    return float_or_none(first_value(mmcif, "_refine.ls_d_res_high")) or float_or_none(
        first_value(mmcif, "_em_3d_reconstruction.resolution")
    )


def entry_methods(mmcif: dict[str, Any]) -> str | None:
    methods = as_list(mmcif.get("_exptl.method"))
    return "; ".join(str(method) for method in methods if method not in {".", "?"}) or None


def r_factor(mmcif: dict[str, Any], key: str) -> float | None:
    return float_or_none(first_value(mmcif, key))


def chem_comp_descriptors(mmcif: dict[str, Any]) -> dict[str, dict[str, str | None]]:
    comp_ids = as_list(mmcif.get("_pdbx_chem_comp_descriptor.comp_id"))
    types = as_list(mmcif.get("_pdbx_chem_comp_descriptor.type"))
    descriptors = as_list(mmcif.get("_pdbx_chem_comp_descriptor.descriptor"))
    result: dict[str, dict[str, str | None]] = {}
    for comp_id, dtype, descriptor in zip(comp_ids, types, descriptors, strict=False):
        if descriptor in {".", "?", None}:
            continue
        comp = result.setdefault(str(comp_id).upper(), {"smiles": None, "inchi_key": None})
        dtype_upper = str(dtype).upper()
        if dtype_upper == "SMILES_CANONICAL":
            comp["smiles"] = str(descriptor)
        elif dtype_upper == "SMILES" and not comp.get("smiles"):
            comp["smiles"] = str(descriptor)
        elif dtype_upper == "INCHIKEY":
            comp["inchi_key"] = str(descriptor)
    return result


def chem_comp_table(mmcif: dict[str, Any]) -> dict[str, dict[str, Any]]:
    ids = as_list(mmcif.get("_chem_comp.id"))
    names = as_list(mmcif.get("_chem_comp.name"))
    formulas = as_list(mmcif.get("_chem_comp.formula"))
    weights = as_list(mmcif.get("_chem_comp.formula_weight"))
    table: dict[str, dict[str, Any]] = {}
    for comp_id, name, formula, weight in zip(ids, names, formulas, weights, strict=False):
        table[str(comp_id).upper()] = {
            "ligand_name": None if name in {".", "?"} else str(name),
            "ligand_formula": None if formula in {".", "?"} else str(formula),
            "ligand_mw": float_or_none(weight),
        }
    return table


def protein_chains(structure: Any) -> list[Any]:
    chains = []
    for chain in structure.get_chains():
        if any(is_aa(residue, standard=True) for residue in chain.get_residues()):
            chains.append(chain)
    return chains


def residue_count(chain: Any) -> int:
    return sum(1 for residue in chain.get_residues() if is_aa(residue, standard=True))

