from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

import pandas as pd
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem

from egfr_dockingforge.common.io import write_table


TRANSFORM_ROWS = [
    ("small_substituent_scan", "[*:1]-[H]>>[*:1]-C", "methyl peripheral heavy atoms"),
    ("halogen_scan", "[*:1]-[H]>>[*:1]-F", "fluorine scan"),
    ("halogen_scan", "[*:1]-[H]>>[*:1]-Cl", "chlorine scan"),
    ("heteroatom_swap", "aromatic_C_to_N", "aromatic CH to N where valid"),
    ("solubilizing_tail_tuning", "[*:1]-[H]>>[*:1]-OC", "methoxy scan"),
    ("conservative_bioisostere", "[*:1]-[H]>>[*:1]-C#N", "nitrile scan"),
    ("brics_recombination", "BRICS_decomposition_seed_fragments", "BRICS-fragment recombination placeholder"),
]


def write_transformation_library(paths: dict[str, Path]) -> pd.DataFrame:
    rows = []
    for i, (klass, rule, desc) in enumerate(TRANSFORM_ROWS, start=1):
        rows.append(
            {
                "transformation_id": f"tr_{i:03d}",
                "transformation_class": klass,
                "smarts_or_rule": rule,
                "description": desc,
                "allowed_attachment_context": "unprotected peripheral atom",
                "disallowed_context": "covalent warhead, metal, invalid valence, protected key interaction atom",
                "medchem_rationale": "conservative analog enumeration for tool-verified Stage 9 screening",
                "source": "rdkit_stage9_transform_library",
                "risk_level": "low" if klass != "brics_recombination" else "medium",
                "enabled_by_default": klass != "brics_recombination",
                "warnings_json": json.dumps([]),
            }
        )
    out = pd.DataFrame(rows)
    write_table(paths["processed"] / "transformation_library.parquet", out)
    write_table(paths["processed"] / "transformation_library.csv", out)
    return out


def _canonical(smiles: str) -> str | None:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    try:
        Chem.SanitizeMol(mol)
    except Exception:
        return None
    return Chem.MolToSmiles(mol, isomericSmiles=True)


def _fp(smiles: str):
    mol = Chem.MolFromSmiles(smiles)
    return AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048) if mol is not None else None


def _attach_to_atom(parent: str, atom_idx: int, fragment: str) -> str | None:
    mol = Chem.RWMol(Chem.AddHs(Chem.MolFromSmiles(parent)))
    atom = mol.GetAtomWithIdx(atom_idx)
    h_neighbor = next((n.GetIdx() for n in atom.GetNeighbors() if n.GetAtomicNum() == 1), None)
    if h_neighbor is None:
        return None
    mol.RemoveAtom(h_neighbor)
    attach_idx = atom_idx if h_neighbor > atom_idx else atom_idx - 1
    frag = Chem.MolFromSmiles(fragment)
    combo = Chem.CombineMols(Chem.RemoveHs(mol), frag)
    rw = Chem.RWMol(combo)
    new_idx = Chem.RemoveHs(mol).GetNumAtoms()
    rw.AddBond(attach_idx, new_idx, Chem.BondType.SINGLE)
    try:
        out = Chem.RemoveHs(rw.GetMol())
        Chem.SanitizeMol(out)
    except Exception:
        return None
    return Chem.MolToSmiles(out, isomericSmiles=True)


def enumerate_rule_based_analogs(
    seeds: pd.DataFrame,
    edit_sites: pd.DataFrame,
    config: dict[str, Any],
    paths: dict[str, Path],
    strategy_name: str = "rdkit_rule_based",
) -> pd.DataFrame:
    max_per_seed = int(config["transforms"]["max_analogs_per_seed_per_strategy"])
    fragments = [
        ("small_substituent_scan", "C"),
        ("halogen_scan", "F"),
        ("halogen_scan", "Cl"),
        ("solubilizing_tail_tuning", "OC"),
        ("conservative_bioisostere", "C#N"),
    ]
    known = set()
    rows = []
    for seed in seeds.to_dict("records"):
        parent_fp = _fp(seed["standard_smiles"])
        n = 0
        sites = edit_sites[(edit_sites["seed_id"] == seed["seed_id"]) & ~edit_sites["protected_region_flag"].astype(bool)]
        for site in sites.to_dict("records"):
            atom_idx = int(site["attachment_atom_idx"])
            for klass, frag in fragments:
                if n >= max_per_seed:
                    break
                analog = _canonical(_attach_to_atom(seed["standard_smiles"], atom_idx, frag) or "")
                if not analog:
                    continue
                key = (seed["seed_id"], analog)
                status = "unique"
                if key in known or analog == seed["standard_smiles"]:
                    status = "duplicate"
                known.add(key)
                fp = _fp(analog)
                tanimoto = DataStructs.TanimotoSimilarity(parent_fp, fp) if parent_fp is not None and fp is not None else 0.0
                digest = hashlib.sha1(f"{seed['seed_id']}|{analog}|{strategy_name}".encode()).hexdigest()[:12]
                rows.append(
                    {
                        "analog_id": f"analog_{digest}",
                        "proposal_id": f"proposal_{digest}",
                        "iteration_id": "iter_001",
                        "strategy_name": strategy_name,
                        "seed_id": seed["seed_id"],
                        "parent_molecule_id": seed["molecule_id"],
                        "parent_smiles": seed["standard_smiles"],
                        "analog_smiles": analog,
                        "standard_smiles": analog,
                        "inchi_key": Chem.MolToInchiKey(Chem.MolFromSmiles(analog)),
                        "transformation_class": klass,
                        "edit_site_id": site["edit_site_id"],
                        "generated_by": "rdkit_deterministic_transform",
                        "source": "stage9_rdkit_rule_based",
                        "uniqueness_status": status,
                        "novelty_status": "analog_of_stage8_seed",
                        "parent_tanimoto": tanimoto,
                        "closest_known_egfr_ligand": seed["molecule_id"],
                        "warnings_json": json.dumps([]),
                    }
                )
                if status == "unique":
                    n += 1
    out = pd.DataFrame(rows).drop_duplicates(["seed_id", "standard_smiles", "strategy_name"])
    write_table(paths["processed"] / "analog_candidates.parquet", out)
    write_table(paths["processed"] / "analog_candidates.csv", out)
    return out
