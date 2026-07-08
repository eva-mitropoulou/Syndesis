from __future__ import annotations

from pathlib import Path
from typing import Any

import pandas as pd
from Bio.PDB import PDBParser

from egfr_dockingforge.stage2.schemas import POCKET_MAPPING_COLUMNS


EGFR_RESIDUE_NAME_HINTS: dict[int, set[str]] = {
    745: {"LYS"},
    762: {"GLU"},
    790: {"THR", "MET"},
    793: {"MET"},
    797: {"CYS"},
    855: {"ASP"},
    856: {"PHE"},
    857: {"GLY"},
}


def parse_receptor(path: str | Path) -> Any:
    parser = PDBParser(QUIET=True)
    return parser.get_structure(Path(path).stem, str(path))


def first_chain(structure: Any) -> Any:
    return next(structure.get_chains())


def residue_by_auth_seq(chain: Any) -> dict[int, Any]:
    return {int(residue.id[1]): residue for residue in chain.get_residues()}


def residue_matches_uniprot_hint(residue: Any | None, uniprot_residue: int) -> bool:
    expected = EGFR_RESIDUE_NAME_HINTS.get(int(uniprot_residue))
    if residue is None:
        return False
    if not expected:
        return True
    return residue.get_resname().strip().upper() in expected


def infer_uniprot_to_auth_offset(residues: dict[int, Any]) -> int:
    best_offset = 0
    best_score = -1
    for offset in range(-80, 81):
        score = 0
        for uniprot_residue, expected_names in EGFR_RESIDUE_NAME_HINTS.items():
            residue = residues.get(uniprot_residue + offset)
            if residue is not None and residue.get_resname().strip().upper() in expected_names:
                score += 1
        if score > best_score:
            best_score = score
            best_offset = offset
    return best_offset if best_score >= 3 else 0


def resolve_uniprot_residue(residues: dict[int, Any], uniprot_residue: int) -> tuple[Any | None, int | None, str]:
    direct = residues.get(int(uniprot_residue))
    if residue_matches_uniprot_hint(direct, int(uniprot_residue)):
        return direct, int(uniprot_residue), "direct_auth_seq_id"
    offset = infer_uniprot_to_auth_offset(residues)
    mapped_auth = int(uniprot_residue) + offset
    mapped = residues.get(mapped_auth)
    if offset and residue_matches_uniprot_hint(mapped, int(uniprot_residue)):
        return mapped, mapped_auth, f"egfr_sequence_offset_{offset:+d}"
    return direct, int(uniprot_residue) if direct is not None else None, "unresolved"


def atom_names(residue: Any | None) -> list[str]:
    if residue is None:
        return []
    return sorted(atom.get_name().strip() for atom in residue.get_atoms())


def map_pocket_residues(row: pd.Series, config: dict[str, Any]) -> list[dict[str, Any]]:
    receptor_id = row["receptor_id"]
    structure = parse_receptor(row["receptor_file_path"])
    chain = first_chain(structure)
    residues = residue_by_auth_seq(chain)
    rows: list[dict[str, Any]] = []
    alignment = set(config["pocket"].get("alignment_residues_uniprot", []))
    for uniprot_residue in config["pocket"]["key_residues_uniprot"]:
        residue, auth_seq_id, mapping_method = resolve_uniprot_residue(residues, int(uniprot_residue))
        names = atom_names(residue)
        warning = ""
        if residue is None:
            warning = "Residue missing after EGFR sequence-offset mapping."
        elif "CA" not in names:
            warning = "Residue lacks CA atom."
        rows.append(
            {
                "receptor_id": receptor_id,
                "pdb_id": row["pdb_id"],
                "auth_asym_id": row["auth_asym_id"],
                "uniprot_residue_number": int(uniprot_residue),
                "auth_seq_id": auth_seq_id,
                "residue_name": residue.get_resname().strip() if residue is not None else None,
                "klifs_position": mapping_method,
                "atom_names_present": ",".join(names),
                "residue_complete_flag": residue is not None and "CA" in names,
                "used_for_alignment_flag": int(uniprot_residue) in alignment and residue is not None and "CA" in names,
                "warning": warning,
            }
        )
    return rows


def mapping_frame(rows: list[dict[str, Any]]) -> pd.DataFrame:
    return pd.DataFrame(rows, columns=POCKET_MAPPING_COLUMNS)
