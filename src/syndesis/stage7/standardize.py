from __future__ import annotations

import hashlib
import json
from typing import Any

from rdkit import Chem
from rdkit.Chem.MolStandardize import rdMolStandardize


def molecule_id_from_inchikey(inchi_key: str, smiles: str) -> str:
    token = inchi_key or hashlib.sha1(smiles.encode("utf-8")).hexdigest()[:16]
    return "mol_" + token.replace("-", "_").lower()


def standardize_smiles(raw_smiles: str | None, raw_record_id: str, source: str) -> dict[str, Any]:
    warnings: list[str] = []
    if not raw_smiles:
        return _failed(raw_record_id, source, raw_smiles, "missing_smiles")
    mol = Chem.MolFromSmiles(str(raw_smiles), sanitize=True)
    if mol is None:
        return _failed(raw_record_id, source, raw_smiles, "invalid_smiles")
    mixture = len(Chem.GetMolFrags(mol)) > 1
    metal = any(atom.GetAtomicNum() in {3, 4, 11, 12, 13, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 40, 46, 47, 78, 79} for atom in mol.GetAtoms())
    parent = rdMolStandardize.FragmentParent(mol)
    normalizer = rdMolStandardize.Normalizer()
    mol = normalizer.normalize(parent)
    uncharger = rdMolStandardize.Uncharger()
    mol = uncharger.uncharge(mol)
    tautomer = rdMolStandardize.TautomerEnumerator().Canonicalize(mol)
    smiles = Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True)
    tautomer_smiles = Chem.MolToSmiles(tautomer, canonical=True, isomericSmiles=True)
    inchi_key = Chem.MolToInchiKey(mol)
    if mixture:
        warnings.append("mixture_or_salt_detected_largest_fragment_kept")
    if metal:
        warnings.append("metal_detected")
    return {
        "raw_record_id": raw_record_id,
        "source": source,
        "raw_smiles": raw_smiles,
        "standard_smiles": smiles,
        "canonical_tautomer_smiles": tautomer_smiles,
        "inchi_key": inchi_key,
        "inchi_key_connectivity": inchi_key.split("-")[0] if inchi_key else "",
        "salt_removed_flag": bool(mixture),
        "mixture_flag": bool(mixture),
        "metal_flag": bool(metal),
        "sanitization_status": "success",
        "standardization_status": "success",
        "warnings_json": json.dumps(warnings),
    }


def _failed(raw_record_id: str, source: str, raw_smiles: str | None, reason: str) -> dict[str, Any]:
    return {
        "raw_record_id": raw_record_id, "source": source, "raw_smiles": raw_smiles or "",
        "standard_smiles": "", "canonical_tautomer_smiles": "", "inchi_key": "",
        "inchi_key_connectivity": "", "salt_removed_flag": False, "mixture_flag": False,
        "metal_flag": False, "sanitization_status": "failed", "standardization_status": "failed",
        "warnings_json": json.dumps([reason]),
    }
