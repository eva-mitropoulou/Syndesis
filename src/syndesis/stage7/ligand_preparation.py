from __future__ import annotations

import json

import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem


def prepare_ligands(master: pd.DataFrame, subsets: pd.DataFrame, config: dict, paths: dict) -> pd.DataFrame:
    rows = []
    out_dir = paths["processed"] / "prepared_ligands"
    selected = subsets[subsets["include_in_subset"]].drop_duplicates("molecule_id")
    source = master.set_index("molecule_id")
    for _, item in selected.iterrows():
        molrow = source.loc[item["molecule_id"]]
        mol = Chem.AddHs(Chem.MolFromSmiles(molrow["standard_smiles"]))
        status = "success"
        warnings = []
        if AllChem.EmbedMolecule(mol, randomSeed=707) != 0:
            status = "failed"
            warnings.append("rdkit_embed_failed")
        else:
            AllChem.UFFOptimizeMolecule(mol, maxIters=200)
        prepared_id = f"prep_{item['molecule_id']}"
        sdf_path = out_dir / f"{prepared_id}.sdf"
        if status == "success":
            writer = Chem.SDWriter(str(sdf_path))
            writer.write(mol)
            writer.close()
        rows.append(
            {
                "prepared_ligand_id": prepared_id, "molecule_id": item["molecule_id"], "source": molrow["source"],
                "screening_subset": item["screening_subset"], "standard_smiles": molrow["standard_smiles"],
                "prepared_smiles": Chem.MolToSmiles(Chem.RemoveHs(mol), canonical=True) if status == "success" else "",
                "protonation_state_id": "stage7_rdkit_default", "tautomer_state_id": "canonical_tautomer",
                "stereoisomer_id": "input_stereo", "conformer_id": "conf_0", "sdf_path": str(sdf_path) if status == "success" else "",
                "mol2_path_if_available": "", "pdbqt_path_if_available": "", "ligand_prep_tool": "rdkit_etkdg_uff",
                "pH_target": config["ligand_preparation"]["ph_target"], "pH_window": config["ligand_preparation"]["ph_window"],
                "charge_model": "rdkit_formal_charge", "num_conformers_generated": int(status == "success"),
                "preparation_status": status, "warnings_json": json.dumps(warnings),
            }
        )
    return pd.DataFrame(rows)
