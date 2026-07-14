from __future__ import annotations

import json

from rdkit import Chem
from rdkit.Chem import Crippen, Descriptors, Lipinski, QED, rdMolDescriptors


WARHEAD_SMARTS = [Chem.MolFromSmarts(x) for x in ["C=CC(=O)N", "C=CC(=O)", "C#CC(=O)", "NS(=O)(=O)F"]]


def flags(smiles: str, source: str, config: dict) -> dict:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return {"invalid_structure_flag": True, "include_in_screening_library": False, "exclusion_reason": "invalid_structure"}
    props = {
        "mw": float(Descriptors.MolWt(mol)),
        "heavy_atom_count": int(mol.GetNumHeavyAtoms()),
        "clogp": float(Crippen.MolLogP(mol)),
        "tpsa": float(rdMolDescriptors.CalcTPSA(mol)),
        "hbd": int(Lipinski.NumHDonors(mol)),
        "hba": int(Lipinski.NumHAcceptors(mol)),
        "rotatable_bonds": int(Lipinski.NumRotatableBonds(mol)),
        "ring_count": int(rdMolDescriptors.CalcNumRings(mol)),
        "aromatic_ring_count": int(rdMolDescriptors.CalcNumAromaticRings(mol)),
        "formal_charge": int(sum(a.GetFormalCharge() for a in mol.GetAtoms())),
        "fraction_sp3": float(rdMolDescriptors.CalcFractionCSP3(mol)),
        "qed": float(QED.qed(mol)),
    }
    metal = any(a.GetAtomicNum() in {3, 4, 11, 12, 13, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 46, 47, 78, 79} for a in mol.GetAtoms())
    warhead = any(mol.HasSubstructMatch(patt) for patt in WARHEAD_SMARTS if patt is not None)
    macro = props["ring_count"] > 0 and any(len(ring) >= 12 for ring in mol.GetRingInfo().AtomRings())
    cfg = config["filters"]
    prop_pass = (
        cfg["mw_min"] <= props["mw"] <= cfg["mw_max"]
        and cfg["heavy_atom_min"] <= props["heavy_atom_count"] <= cfg["heavy_atom_max"]
        and cfg["clogp_min"] <= props["clogp"] <= cfg["clogp_max"]
        and cfg["tpsa_min"] <= props["tpsa"] <= cfg["tpsa_max"]
        and props["hbd"] <= cfg["hbd_max"]
        and props["hba"] <= cfg["hba_max"]
        and props["rotatable_bonds"] <= cfg["rotatable_bonds_max"]
        and cfg["formal_charge_min"] <= props["formal_charge"] <= cfg["formal_charge_max"]
    )
    hard_scope = prop_pass and not metal and not warhead and not (macro and cfg["exclude_macrocycles"])
    risks = []
    if props["clogp"] > 5:
        risks.append("high_lipophilicity")
    if props["qed"] < 0.25:
        risks.append("low_qed")
    return {
        **props,
        "pains_flag": False,
        "brenk_flag": False,
        "reactive_flag": bool(warhead),
        "aggregator_risk_flag": props["clogp"] > 5.5,
        "covalent_warhead_flag": bool(warhead),
        "egfr_cys797_warhead_flag": bool(warhead),
        "allosteric_scope_flag": False,
        "macrocycle_flag": bool(macro),
        "metal_flag": bool(metal),
        "mixture_flag": False,
        "invalid_structure_flag": False,
        "property_window_pass": bool(prop_pass),
        "hard_scope_pass": bool(hard_scope),
        "soft_medchem_pass": len(risks) == 0,
        "include_in_screening_library": bool(hard_scope and source not in {"chembl_known_ligand", "bindingdb_known_ligand", "stage1_native_ligand"}),
        "exclusion_reason": "" if hard_scope else "hard_scope_or_property_window_failure",
        "risk_flags_json": json.dumps(risks),
    }
