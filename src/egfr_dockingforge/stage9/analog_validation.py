from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pandas as pd
from rdkit import Chem
from rdkit.Chem import Crippen, Descriptors, QED, rdMolDescriptors

from egfr_dockingforge.common.io import write_table

WARHEAD_SMARTS = [
    Chem.MolFromSmarts("C=CC(=O)N"),
    Chem.MolFromSmarts("C=CC(=O)O"),
    Chem.MolFromSmarts("C#CC(=O)N"),
    Chem.MolFromSmarts("S(=O)(=O)F"),
]
REACTIVE_SMARTS = [
    Chem.MolFromSmarts("[N+](=O)[O-]"),
    Chem.MolFromSmarts("[CX3](=O)Cl"),
    Chem.MolFromSmarts("[SH]"),
]


def _has_match(mol: Chem.Mol, patterns: list[Chem.Mol | None]) -> bool:
    return any(patt is not None and mol.HasSubstructMatch(patt) for patt in patterns)


def validate_analog_batch(candidates: pd.DataFrame, config: dict[str, Any], paths: dict[str, Path]) -> pd.DataFrame:
    rows = []
    for row in candidates.to_dict("records"):
        mol = Chem.MolFromSmiles(row.get("standard_smiles", ""))
        if mol is None:
            rows.append(
                {
                    "analog_id": row["analog_id"],
                    "valid_molecule_flag": False,
                    "standardization_status": "invalid_smiles",
                    "hard_scope_pass": False,
                    "covalent_warhead_flag": False,
                    "reactive_flag": False,
                    "pains_flag": False,
                    "brenk_flag": False,
                    "property_window_pass": False,
                    "mw": None,
                    "clogp": None,
                    "tpsa": None,
                    "hbd": None,
                    "hba": None,
                    "rotatable_bonds": None,
                    "qed": None,
                    "sa_score_if_available": None,
                    "medchem_risk_score": 1.0,
                    "validation_status": "rejected",
                    "rejection_reason": "invalid_smiles",
                    "warnings_json": json.dumps([]),
                }
            )
            continue
        mw = Descriptors.MolWt(mol)
        clogp = Crippen.MolLogP(mol)
        tpsa = rdMolDescriptors.CalcTPSA(mol)
        hbd = rdMolDescriptors.CalcNumHBD(mol)
        hba = rdMolDescriptors.CalcNumHBA(mol)
        rot = rdMolDescriptors.CalcNumRotatableBonds(mol)
        qed = QED.qed(mol)
        warhead = _has_match(mol, WARHEAD_SMARTS)
        reactive = _has_match(mol, REACTIVE_SMARTS)
        property_pass = 180 <= mw <= 650 and -1.0 <= clogp <= 6.5 and tpsa <= 160 and rot <= 12
        hard_scope = not warhead and not reactive and property_pass and not any(a.GetAtomicNum() > 53 for a in mol.GetAtoms())
        risk = min(1.0, (0.3 if warhead else 0) + (0.3 if reactive else 0) + (0 if property_pass else 0.3) + max(clogp - 5, 0) * 0.1)
        reason = "" if hard_scope else ("covalent_warhead" if warhead else "reactive_or_property_scope_failure")
        rows.append(
            {
                "analog_id": row["analog_id"],
                "valid_molecule_flag": True,
                "standardization_status": "standardized",
                "hard_scope_pass": bool(hard_scope),
                "covalent_warhead_flag": bool(warhead),
                "reactive_flag": bool(reactive),
                "pains_flag": False,
                "brenk_flag": False,
                "property_window_pass": bool(property_pass),
                "mw": mw,
                "clogp": clogp,
                "tpsa": tpsa,
                "hbd": hbd,
                "hba": hba,
                "rotatable_bonds": rot,
                "qed": qed,
                "sa_score_if_available": None,
                "medchem_risk_score": risk,
                "validation_status": "accepted_for_screening" if hard_scope else "rejected",
                "rejection_reason": reason,
                "warnings_json": json.dumps([]),
            }
        )
    out = pd.DataFrame(rows)
    write_table(paths["processed"] / "analog_validation.parquet", out)
    write_table(paths["processed"] / "analog_validation.csv", out)
    return out
