from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from syndesis.stage2.pocket_mapping import first_chain, parse_receptor, residue_by_auth_seq, resolve_uniprot_residue
from syndesis.stage2.state_labels import normalize_state_label
from syndesis.stage2.schemas import FEATURE_COLUMNS


def atom_coord(residue: Any | None, atom_name: str) -> np.ndarray | None:
    if residue is None:
        return None
    for atom in residue.get_atoms():
        if atom.get_name().strip() == atom_name:
            return np.asarray(atom.coord, dtype=float)
    return None


def min_atom_distance(residue_a: Any | None, atom_a: str, residue_b: Any | None, atom_prefixes_b: tuple[str, ...]) -> float | None:
    coord_a = atom_coord(residue_a, atom_a)
    if coord_a is None or residue_b is None:
        return None
    distances = []
    for atom in residue_b.get_atoms():
        if atom.get_name().strip().startswith(atom_prefixes_b):
            distances.append(float(np.linalg.norm(coord_a - np.asarray(atom.coord, dtype=float))))
    return min(distances) if distances else None


def pdb_heavy_atom_coords(path: str | Path) -> np.ndarray:
    coords = []
    with Path(path).open("r", encoding="utf-8") as handle:
        for line in handle:
            if not line.startswith(("ATOM", "HETATM")):
                continue
            element = line[76:78].strip() or line[12:16].strip()[0]
            if element.upper() == "H":
                continue
            coords.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
    return np.asarray(coords, dtype=float)


def count_pdb_records(path: str | Path, prefixes: tuple[str, ...]) -> int:
    with Path(path).open("r", encoding="utf-8") as handle:
        return sum(1 for line in handle if line.startswith(prefixes))


def centroid_and_rg(coords: np.ndarray) -> tuple[np.ndarray | None, float | None]:
    if coords.size == 0:
        return None, None
    centroid = coords.mean(axis=0)
    rg = float(np.sqrt(((coords - centroid) ** 2).sum(axis=1).mean()))
    return centroid, rg


def residue_complete(residues: dict[int, Any], residue_number: int) -> bool:
    residue, _, _ = resolve_uniprot_residue(residues, residue_number)
    return residue is not None and atom_coord(residue, "CA") is not None


def preselection_score(row: pd.Series, completeness: float, warnings: list[str]) -> float:
    score = float(row.get("quality_score") or 0.0)
    score += completeness * 10.0
    if row.get("quality_tier") == "Tier A":
        score += 5.0
    if row.get("quality_tier") == "Tier C":
        score -= 10.0
    score -= len(warnings) * 1.0
    return round(max(0.0, score), 3)


def feature_record(row: pd.Series, config: dict[str, Any]) -> dict[str, Any]:
    receptor_id = row["receptor_id"]
    warnings: list[str] = []
    structure = parse_receptor(row["receptor_file_path"])
    chain = first_chain(structure)
    residues = residue_by_auth_seq(chain)
    lys, lys_auth, _ = resolve_uniprot_residue(residues, 745)
    glu, glu_auth, _ = resolve_uniprot_residue(residues, 762)
    phe, _, _ = resolve_uniprot_residue(residues, 856)
    phe_ca = atom_coord(phe, "CA")
    saltbridge = min_atom_distance(lys, "NZ", glu, ("OE",))
    if saltbridge is None:
        warnings.append("Lys745-Glu762 salt-bridge distance unavailable after EGFR sequence-offset mapping.")
    if phe_ca is None:
        warnings.append("Phe856 CA position unavailable under fallback mapping.")

    ligand_coords = pdb_heavy_atom_coords(row["native_ligand_sdf_path"])
    centroid, rg = centroid_and_rg(ligand_coords)
    if centroid is None:
        warnings.append("Native ligand heavy-atom centroid unavailable.")
        centroid = np.array([np.nan, np.nan, np.nan])
    state = normalize_state_label(row)
    if state["state_stratum"] == "unknown_state":
        warnings.append("KinCore/KLIFS state labels unavailable; state_stratum set to unknown_state.")
    completeness = float(row.get("active_site_completeness_score") or 0.0)
    record = {
        "receptor_id": receptor_id,
        "complex_id": row["complex_id"],
        "pdb_id": row["pdb_id"],
        "auth_asym_id": row["auth_asym_id"],
        "ligand_comp_id": row["ligand_comp_id"],
        "quality_tier": row["quality_tier"],
        "quality_score": row["quality_score"],
        "resolution_angstrom": row["resolution_angstrom"],
        "mutation_flag": bool(row.get("mutation_flag")) if pd.notna(row.get("mutation_flag")) else False,
        "mutation_list": "" if pd.isna(row.get("mutation_list")) else row.get("mutation_list"),
        "kincore_activity_label": state["kincore_activity_label"],
        "dfg_state": row.get("dfg_state"),
        "dfg_dihedral_cluster": state["dfg_dihedral_cluster"],
        "chelix_state": row.get("chelix_state"),
        "saltbridge_state": "saltbridge_in" if saltbridge is not None and saltbridge <= 4.0 else None,
        "hrd_state": row.get("hrd_state"),
        "activation_loop_state": row.get("activation_loop_state"),
        "klifs_structure_id": row.get("klifs_structure_id"),
        "klifs_pocket_id": row.get("klifs_pocket_id"),
        "ligand_class_if_known": row.get("inhibitor_type_if_known"),
        "pocket_ca_rmsd_to_reference": None,
        "hinge_region_rmsd": None,
        "dfg_region_rmsd": None,
        "c_helix_proxy_distance": saltbridge,
        "lys745_glu762_nz_oe_min_distance": saltbridge,
        "dfg_phe856_position_x": float(phe_ca[0]) if phe_ca is not None else None,
        "dfg_phe856_position_y": float(phe_ca[1]) if phe_ca is not None else None,
        "dfg_phe856_position_z": float(phe_ca[2]) if phe_ca is not None else None,
        "gatekeeper_thr790_complete": residue_complete(residues, 790),
        "hinge_met793_complete": residue_complete(residues, 793),
        "asp855_complete": residue_complete(residues, 855),
        "phe856_complete": residue_complete(residues, 856),
        "pocket_volume_if_available": None,
        "native_ligand_centroid_x": float(centroid[0]),
        "native_ligand_centroid_y": float(centroid[1]),
        "native_ligand_centroid_z": float(centroid[2]),
        "native_ligand_heavy_atom_count": int(len(ligand_coords)),
        "native_ligand_radius_of_gyration": rg,
        "pocket_water_count_5a": count_pdb_records(row.get("pocket_waters_path", ""), ("HETATM",)) if isinstance(row.get("pocket_waters_path"), str) and Path(row.get("pocket_waters_path")).exists() else 0,
        "active_site_completeness_score": completeness,
        "receptor_file_path": row["receptor_clean_path"],
        "native_complex_path": row["native_complex_path"],
        "native_ligand_sdf_path": row["native_ligand_sdf_path"],
        "state_stratum": state["state_stratum"],
        "warnings_json": None,
    }
    record["receptor_preselection_score"] = preselection_score(row, completeness, warnings)
    record["warnings_json"] = json.dumps(warnings, sort_keys=True)
    return {column: record.get(column) for column in FEATURE_COLUMNS}


def features_frame(records: list[dict[str, Any]]) -> pd.DataFrame:
    return pd.DataFrame(records, columns=FEATURE_COLUMNS)
