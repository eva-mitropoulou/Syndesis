from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pandas as pd
from rdkit import Chem

from egfr_dockingforge.common.io import write_table


def detect_edit_sites(seeds: pd.DataFrame, config: dict[str, Any], paths: dict[str, Path]) -> pd.DataFrame:
    rows = []
    max_sites = int(config["edit_sites"]["max_sites_per_seed"])
    allowed = json.dumps(config["transforms"]["enabled_classes"])
    disallowed = json.dumps(["covalent_warhead_introduction", "core_demolition", "allosteric_design"])
    for seed in seeds.to_dict("records"):
        mol = Chem.MolFromSmiles(seed["standard_smiles"])
        if mol is None:
            rows.append(
                {
                    "edit_site_id": f"{seed['seed_id']}_no_site",
                    "seed_id": seed["seed_id"],
                    "molecule_id": seed["molecule_id"],
                    "atom_indices_json": "[]",
                    "attachment_atom_idx": -1,
                    "editable_region_type": "none",
                    "protected_region_flag": True,
                    "protected_reason": "invalid_seed_smiles",
                    "solvent_exposed_proxy": 0.0,
                    "interaction_region": "unknown",
                    "allowed_transformation_classes_json": "[]",
                    "disallowed_transformation_classes_json": disallowed,
                    "evidence_json": "{}",
                    "warnings_json": json.dumps(["missing_edit_sites"]),
                }
            )
            continue
        candidates = []
        for atom in mol.GetAtoms():
            if atom.GetAtomicNum() in {6, 7, 8, 16} and atom.GetDegree() <= 3:
                candidates.append(atom.GetIdx())
        for j, idx in enumerate(candidates[:max_sites], start=1):
            atom = mol.GetAtomWithIdx(idx)
            protected = atom.GetIsAromatic() and atom.GetDegree() >= 3
            rows.append(
                {
                    "edit_site_id": f"{seed['seed_id']}_site_{j:02d}",
                    "seed_id": seed["seed_id"],
                    "molecule_id": seed["molecule_id"],
                    "atom_indices_json": json.dumps([idx]),
                    "attachment_atom_idx": idx,
                    "editable_region_type": "r_group_or_peripheral_atom",
                    "protected_region_flag": bool(protected),
                    "protected_reason": "aromatic_core_branchpoint" if protected else "",
                    "solvent_exposed_proxy": round(1.0 / max(atom.GetDegree(), 1), 3),
                    "interaction_region": "peripheral_or_solvent_exposed_proxy",
                    "allowed_transformation_classes_json": "[]" if protected else allowed,
                    "disallowed_transformation_classes_json": disallowed,
                    "evidence_json": json.dumps({"atom_symbol": atom.GetSymbol(), "degree": atom.GetDegree()}),
                    "warnings_json": json.dumps([]),
                }
            )
    out = pd.DataFrame(rows)
    write_table(paths["processed"] / "edit_sites.parquet", out)
    write_table(paths["processed"] / "edit_sites.csv", out)
    return out
