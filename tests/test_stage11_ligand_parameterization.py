from pathlib import Path

import pandas as pd
import pytest

from egfr_dockingforge.stage11.ligand_parameterization import parameterize_ligands


def test_missing_cgenff_str_fails_explicitly_when_cgenff_backend_selected(tmp_path):
    candidates = pd.DataFrame([{"md_candidate_id":"m1","molecule_id":"mol_a","standard_smiles":"CC","ligand_file":"lig.sdf"}])
    config = {"forcefield":{"ligand_parameterization_backend":"charmm_cgenff_paramchem","cgenff_executable":"cgenff","charmm2gmx_executable":"cgenff_charmm2gmx"}}
    out = parameterize_ligands(candidates, config, {"processed": tmp_path, "user_cgenff_str_dir": tmp_path / "str"})
    assert out["parameterization_status"].item() == "failed_parameterization"
    assert out["rejection_reason"].item() == "missing_required_cgenff_str_file"


def test_gaff_backend_rejects_charmm_protein_stack(tmp_path):
    candidates = pd.DataFrame([{"md_candidate_id":"m1","molecule_id":"mol_a","standard_smiles":"CC","ligand_file":"lig.sdf"}])
    config = {"forcefield":{"ligand_parameterization_backend":"amber_gaff2_acpype","protein_ff":"CHARMM36m"}}
    with pytest.raises(ValueError):
        parameterize_ligands(candidates, config, {"processed": tmp_path, "md_root": tmp_path, "user_cgenff_str_dir": tmp_path / "str"})
