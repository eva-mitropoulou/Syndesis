import pandas as pd
from rdkit import Chem

from egfr_dockingforge.stage9.rdkit_transformations import enumerate_rule_based_analogs, write_transformation_library


def test_rdkit_transformations_make_valid_nonduplicate_molecules(tmp_path):
    seeds = pd.DataFrame([{"seed_id": "seed_001", "molecule_id": "mol1", "standard_smiles": "c1ccccc1"}])
    sites = pd.DataFrame(
        [
            {
                "seed_id": "seed_001",
                "edit_site_id": "seed_001_site_01",
                "attachment_atom_idx": 0,
                "protected_region_flag": False,
            }
        ]
    )
    config = {"transforms": {"max_analogs_per_seed_per_strategy": 3}}
    write_transformation_library({"processed": tmp_path})
    out = enumerate_rule_based_analogs(seeds, sites, config, {"processed": tmp_path})
    assert not out.empty
    assert out["parent_molecule_id"].eq("mol1").all()
    assert all(Chem.MolFromSmiles(smi) is not None for smi in out["standard_smiles"])
