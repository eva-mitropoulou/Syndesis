from __future__ import annotations

import hashlib
import os
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import pandas as pd
import requests
from Bio.PDB import MMCIFParser, NeighborSearch, PDBIO, Select
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from Bio.PDB.Polypeptide import is_aa
from joblib import Parallel, delayed

from syndesis.common.io import ensure_dir, load_yaml, project_root, resolve_path, write_json, write_table

RCSB_SEARCH_URL = "https://search.rcsb.org/rcsbsearch/v2/query"
RCSB_MMCIF_URL = "https://files.rcsb.org/download/{pdb_id}.cif"
RCSB_CCD_URL = "https://files.rcsb.org/ligands/view/{comp_id}.cif"
WATER_RESNAMES = {"HOH", "WAT", "DOD", "H2O"}
COMMON_NON_LIGANDS = {
    "ACT",
    "BME",
    "CL",
    "EDO",
    "FMT",
    "GOL",
    "MG",
    "NA",
    "NO3",
    "PEG",
    "PO4",
    "SO4",
}


@dataclass(frozen=True)
class LigandHit:
    chain_id: str
    residue_id: tuple[str, int, str]
    resname: str
    heavy_atoms: int
    carbon_atoms: int
    min_protein_distance: float | None


class ProteinSelect(Select):
    def __init__(self, chain_id: str | None = None) -> None:
        self.chain_id = chain_id

    def accept_chain(self, chain: Any) -> bool:
        return self.chain_id is None or chain.id == self.chain_id

    def accept_residue(self, residue: Any) -> bool:
        return bool(is_aa(residue, standard=True))


class LigandSelect(Select):
    def __init__(self, chain_id: str, residue_id: tuple[str, int, str]) -> None:
        self.chain_id = chain_id
        self.residue_id = residue_id

    def accept_chain(self, chain: Any) -> bool:
        return chain.id == self.chain_id

    def accept_residue(self, residue: Any) -> bool:
        return residue.id == self.residue_id


def normalize_pdb_id(value: str) -> str:
    cleaned = "".join(char for char in str(value).strip().upper() if char.isalnum())
    if len(cleaned) != 4:
        raise ValueError(f"Invalid PDB ID {value!r}")
    return cleaned


def search_payload(stage_config: dict[str, Any]) -> dict[str, Any]:
    query = stage_config["rcsb_query"]
    nodes: list[dict[str, Any]] = [
        {
            "type": "terminal",
            "service": "text",
            "parameters": {
                "attribute": "rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers.database_accession",
                "operator": "exact_match",
                "value": query.get("uniprot_accession", "P00533"),
            },
        },
        {
            "type": "terminal",
            "service": "text",
            "parameters": {
                "attribute": "exptl.method",
                "operator": "exact_match",
                "value": query.get("experimental_method", "X-RAY DIFFRACTION"),
            },
        },
        {
            "type": "terminal",
            "service": "text",
            "parameters": {
                "attribute": "rcsb_entry_info.resolution_combined",
                "operator": "less_or_equal",
                "value": float(query.get("max_resolution_angstrom", 2.8)),
            },
        },
        {
            "type": "terminal",
            "service": "text",
            "parameters": {
                "attribute": "rcsb_entry_info.nonpolymer_entity_count",
                "operator": "greater",
                "value": 0,
            },
        },
    ]
    return {
        "query": {"type": "group", "logical_operator": "and", "nodes": nodes},
        "return_type": "entry",
        "request_options": {
            "paginate": {"start": 0, "rows": int(query.get("search_rows", 250))},
            "sort": [{"sort_by": "rcsb_entry_info.resolution_combined", "direction": "asc"}],
            "scoring_strategy": "combined",
        },
    }


def query_rcsb(stage_config: dict[str, Any]) -> list[str]:
    manual = stage_config.get("manual_pdb_ids") or []
    if manual:
        return list(dict.fromkeys(normalize_pdb_id(pdb_id) for pdb_id in manual))

    response = requests.post(RCSB_SEARCH_URL, json=search_payload(stage_config), timeout=90)
    response.raise_for_status()
    payload = response.json()
    ids = [normalize_pdb_id(row["identifier"]) for row in payload.get("result_set", [])]
    return list(dict.fromkeys(ids))[: int(stage_config["rcsb_query"].get("max_structures", 40))]


def download_mmcif(pdb_id: str, raw_dir: Path, force: bool) -> Path:
    ensure_dir(raw_dir)
    target = raw_dir / f"{pdb_id.lower()}.cif"
    if target.exists() and target.stat().st_size > 0 and not force:
        return target
    response = requests.get(RCSB_MMCIF_URL.format(pdb_id=pdb_id), timeout=90)
    response.raise_for_status()
    target.write_bytes(response.content)
    return target


def first_value(mmcif: dict[str, Any], key: str) -> Any | None:
    value = mmcif.get(key)
    if value is None:
        return None
    if isinstance(value, list):
        return value[0] if value else None
    return value


def mmcif_resolution(mmcif: dict[str, Any]) -> float | None:
    value = first_value(mmcif, "_refine.ls_d_res_high") or first_value(
        mmcif, "_em_3d_reconstruction.resolution"
    )
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def chem_comp_smiles(mmcif: dict[str, Any]) -> dict[str, str]:
    comp_ids = mmcif.get("_pdbx_chem_comp_descriptor.comp_id", [])
    types = mmcif.get("_pdbx_chem_comp_descriptor.type", [])
    descriptors = mmcif.get("_pdbx_chem_comp_descriptor.descriptor", [])
    if isinstance(comp_ids, str):
        comp_ids = [comp_ids]
        types = [types]
        descriptors = [descriptors]
    smiles: dict[str, str] = {}
    preferred = {"SMILES_CANONICAL", "SMILES"}
    for comp_id, descriptor_type, descriptor in zip(comp_ids, types, descriptors, strict=False):
        if str(descriptor_type).upper() in preferred and descriptor not in {".", "?"}:
            smiles.setdefault(str(comp_id).upper(), str(descriptor))
    return smiles


def chem_comp_identifiers(mmcif: dict[str, Any]) -> dict[str, str | None]:
    types = mmcif.get("_pdbx_chem_comp_descriptor.type", [])
    descriptors = mmcif.get("_pdbx_chem_comp_descriptor.descriptor", [])
    if isinstance(types, str):
        types = [types]
        descriptors = [descriptors]

    smiles: str | None = None
    inchi_key: str | None = None
    for descriptor_type, descriptor in zip(types, descriptors, strict=False):
        dtype = str(descriptor_type).upper()
        if descriptor in {".", "?"}:
            continue
        if dtype == "SMILES_CANONICAL":
            smiles = str(descriptor)
        elif dtype == "SMILES" and smiles is None:
            smiles = str(descriptor)
        elif dtype == "INCHIKEY":
            inchi_key = str(descriptor)
    return {"smiles": smiles, "inchi_key": inchi_key}


def fetch_ccd_identifiers(comp_id: str, ccd_dir: Path) -> dict[str, str | None]:
    ensure_dir(ccd_dir)
    comp_id = comp_id.upper()
    path = ccd_dir / f"{comp_id}.cif"
    if not path.exists() or path.stat().st_size == 0:
        response = requests.get(RCSB_CCD_URL.format(comp_id=comp_id), timeout=60)
        response.raise_for_status()
        path.write_bytes(response.content)
    return chem_comp_identifiers(MMCIF2Dict(str(path)))


def parse_structure(path: Path) -> Any:
    parser = MMCIFParser(QUIET=True)
    return parser.get_structure(path.stem.upper(), str(path))


def residue_atom_counts(residue: Any) -> tuple[int, int]:
    heavy = 0
    carbon = 0
    for atom in residue.get_atoms():
        element = (atom.element or atom.get_name()[0]).upper()
        if element != "H":
            heavy += 1
        if element == "C":
            carbon += 1
    return heavy, carbon


def min_distance_to_protein(residue: Any, protein_atoms: list[Any]) -> float | None:
    ligand_atoms = [atom for atom in residue.get_atoms() if (atom.element or "").upper() != "H"]
    if not ligand_atoms or not protein_atoms:
        return None
    search = NeighborSearch(protein_atoms)
    distances: list[float] = []
    for atom in ligand_atoms:
        neighbors = search.search(atom.coord, 8.0, level="A")
        distances.extend(atom - neighbor for neighbor in neighbors)
    return min(distances) if distances else None


def ligand_hits(structure: Any, stage_config: dict[str, Any]) -> list[LigandHit]:
    curation = stage_config.get("curation", {})
    min_heavy = int(curation.get("min_ligand_heavy_atoms", 12))
    max_distance = float(curation.get("max_ligand_protein_distance_angstrom", 5.0))
    excluded = {name.upper() for name in curation.get("excluded_ligand_resnames", [])}
    excluded.update(COMMON_NON_LIGANDS)
    protein_atoms = [
        atom
        for residue in structure.get_residues()
        if is_aa(residue, standard=True)
        for atom in residue.get_atoms()
        if (atom.element or "").upper() != "H"
    ]

    hits: list[LigandHit] = []
    for residue in structure.get_residues():
        resname = residue.get_resname().strip().upper()
        hetero_flag = str(residue.id[0]).strip()
        if not hetero_flag or resname in WATER_RESNAMES or resname in excluded:
            continue
        if is_aa(residue, standard=True):
            continue
        heavy, carbon = residue_atom_counts(residue)
        if heavy < min_heavy or carbon == 0:
            continue
        distance = min_distance_to_protein(residue, protein_atoms)
        if distance is None or distance > max_distance:
            continue
        hits.append(
            LigandHit(
                chain_id=residue.get_parent().id,
                residue_id=residue.id,
                resname=resname,
                heavy_atoms=heavy,
                carbon_atoms=carbon,
                min_protein_distance=distance,
            )
        )
    hits.sort(key=lambda hit: (hit.heavy_atoms, -float(hit.min_protein_distance or 99.0)), reverse=True)
    return hits[: int(curation.get("max_ligands_per_structure", 1))]


def write_selected_structure_files(
    structure: Any,
    hit: LigandHit,
    structure_dir: Path,
    force: bool,
) -> tuple[Path, Path]:
    ensure_dir(structure_dir)
    receptor_path = structure_dir / "receptor_clean.pdb"
    ligand_path = structure_dir / "native_ligand.pdb"
    io = PDBIO()
    io.set_structure(structure)
    if force or not receptor_path.exists():
        io.save(str(receptor_path), ProteinSelect())
    if force or not ligand_path.exists():
        io.save(str(ligand_path), LigandSelect(hit.chain_id, hit.residue_id))
    return receptor_path, ligand_path


def stable_id(*parts: str) -> str:
    joined = "|".join(parts)
    return hashlib.sha1(joined.encode("utf-8")).hexdigest()[:12]


def process_one_pdb(
    pdb_id: str,
    raw_dir: Path,
    ccd_dir: Path,
    prepared_dir: Path,
    stage_config: dict[str, Any],
    force: bool,
) -> dict[str, Any]:
    result: dict[str, Any] = {
        "pdb_id": pdb_id,
        "status": "failed",
        "structures": [],
        "ligands": [],
        "warnings": [],
    }
    try:
        mmcif_path = download_mmcif(pdb_id, raw_dir, force=force)
        mmcif = MMCIF2Dict(str(mmcif_path))
        structure = parse_structure(mmcif_path)
        smiles_by_comp = chem_comp_smiles(mmcif)
        resolution = mmcif_resolution(mmcif)
        max_resolution = stage_config.get("curation", {}).get("max_resolution_angstrom")
        if max_resolution is not None and resolution is not None and resolution > float(max_resolution):
            result["status"] = "excluded"
            result["warnings"].append(
                f"Resolution {resolution:.2f} A is above cutoff {float(max_resolution):.2f} A."
            )
            return result
        hits = ligand_hits(structure, stage_config)
        if not hits:
            result["status"] = "excluded"
            result["warnings"].append("No qualifying organic protein-bound ligand found.")
            return result

        for hit in hits:
            structure_id = f"{pdb_id}_{hit.chain_id}_{hit.resname}"
            ligand_id = f"{pdb_id}_{hit.chain_id}_{hit.resname}_{hit.residue_id[1]}"
            structure_dir = prepared_dir / structure_id.lower()
            receptor_path, ligand_path = write_selected_structure_files(
                structure=structure,
                hit=hit,
                structure_dir=structure_dir,
                force=force,
            )
            smiles = smiles_by_comp.get(hit.resname)
            inchi_key = None
            if not smiles:
                try:
                    identifiers = fetch_ccd_identifiers(hit.resname, ccd_dir)
                    smiles = identifiers["smiles"]
                    inchi_key = identifiers["inchi_key"]
                except Exception as exc:
                    result["warnings"].append(f"CCD lookup failed for {hit.resname}: {exc}")
            result["structures"].append(
                {
                    "structure_id": structure_id,
                    "pdb_id": pdb_id,
                    "chain_id": hit.chain_id,
                    "mutation_status": "not_curated",
                    "resolution": resolution,
                    "ligand_id": ligand_id,
                    "receptor_state": "pending_stage2",
                    "file_path": str(receptor_path),
                    "include_flag": True,
                    "exclusion_reason": "",
                    "source_mmcif": str(mmcif_path),
                    "native_ligand_file": str(ligand_path),
                    "ligand_resname": hit.resname,
                    "ligand_min_protein_distance": hit.min_protein_distance,
                    "ligand_heavy_atoms": hit.heavy_atoms,
                }
            )
            result["ligands"].append(
                {
                    "ligand_id": ligand_id,
                    "source": "native_cocrystal",
                    "smiles": smiles,
                    "standard_smiles": smiles,
                    "inchi_key": inchi_key,
                    "parent_ligand_id": None,
                    "activity_value": None,
                    "activity_type": None,
                    "activity_units": None,
                    "activity_source": None,
                    "include_flag": True,
                    "pdb_id": pdb_id,
                    "resname": hit.resname,
                    "chain_id": hit.chain_id,
                    "residue_number": hit.residue_id[1],
                    "native_ligand_file": str(ligand_path),
                }
            )
        result["status"] = "complete"
        return result
    except Exception as exc:
        result["warnings"].append(str(exc))
        return result


def empty_structures_frame() -> pd.DataFrame:
    return pd.DataFrame(
        columns=[
            "structure_id",
            "pdb_id",
            "chain_id",
            "mutation_status",
            "resolution",
            "ligand_id",
            "receptor_state",
            "file_path",
            "include_flag",
            "exclusion_reason",
        ]
    )


def empty_ligands_frame() -> pd.DataFrame:
    return pd.DataFrame(
        columns=[
            "ligand_id",
            "source",
            "smiles",
            "standard_smiles",
            "inchi_key",
            "parent_ligand_id",
            "activity_value",
            "activity_type",
            "activity_units",
            "activity_source",
            "include_flag",
        ]
    )


def run_stage1(
    project_config: Path,
    stage_config: Path,
    out_dir: Path,
    workers: int | None = None,
    force: bool = False,
) -> dict[str, Any]:
    root = project_root()
    project = load_yaml(resolve_path(project_config, root))
    stage = load_yaml(resolve_path(stage_config, root))
    out_dir = resolve_path(out_dir, root)
    raw_dir = resolve_path(project["paths"]["raw_structures"], root)
    ccd_dir = resolve_path(project["paths"].get("chemical_components", "data/external/chemical_components"), root)
    prepared_dir = out_dir / "prepared_complexes"
    ensure_dir(raw_dir)
    ensure_dir(ccd_dir)
    ensure_dir(prepared_dir)

    started = time.time()
    pdb_ids = query_rcsb(stage)
    if not pdb_ids:
        raise RuntimeError("No PDB IDs were returned by the Stage 1 configuration.")

    n_workers = workers or int(stage.get("runtime", {}).get("workers", min(8, os.cpu_count() or 1)))
    results = Parallel(n_jobs=max(1, n_workers), prefer="threads")(
        delayed(process_one_pdb)(pdb_id, raw_dir, ccd_dir, prepared_dir, stage, force) for pdb_id in pdb_ids
    )
    structure_rows = [row for result in results for row in result["structures"]]
    ligand_rows = [row for result in results for row in result["ligands"]]
    structures = pd.DataFrame(structure_rows) if structure_rows else empty_structures_frame()
    ligands = pd.DataFrame(ligand_rows) if ligand_rows else empty_ligands_frame()

    write_table(out_dir / "structures.csv", structures)
    write_table(out_dir / "ligands.csv", ligands)
    manifest = {
        "project": project.get("project", {}).get("name", "Syndesis"),
        "stage": "stage1_egfr_cocrystal_benchmark",
        "status": "csv_complete_parquet_pending" if not structures.empty else "no_included_structures",
        "pdb_ids_queried": pdb_ids,
        "included_structures": int(len(structures)),
        "included_ligands": int(len(ligands)),
        "workers": n_workers,
        "elapsed_seconds": round(time.time() - started, 2),
        "structures_table": str(out_dir / "structures.parquet"),
        "ligands_table": str(out_dir / "ligands.parquet"),
        "structures_csv": str(out_dir / "structures.csv"),
        "ligands_csv": str(out_dir / "ligands.csv"),
        "excluded_or_failed": [
            {"pdb_id": result["pdb_id"], "status": result["status"], "warnings": result["warnings"]}
            for result in results
            if result["status"] != "complete"
        ],
    }
    write_json(out_dir / "stage1_manifest.json", manifest)

    structures_table = out_dir / "structures.parquet"
    ligands_table = out_dir / "ligands.parquet"
    try:
        write_table(structures_table, structures)
        write_table(ligands_table, ligands)
    except RuntimeError as exc:
        manifest["status"] = "csv_complete_parquet_failed"
        manifest["parquet_error"] = str(exc)
        manifest["elapsed_seconds"] = round(time.time() - started, 2)
        write_json(out_dir / "stage1_manifest.json", manifest)
        raise
    manifest["status"] = "complete" if not structures.empty else "no_included_structures"
    manifest["elapsed_seconds"] = round(time.time() - started, 2)
    write_json(out_dir / "stage1_manifest.json", manifest)
    return manifest
