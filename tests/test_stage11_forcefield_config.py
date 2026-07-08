from egfr_dockingforge.stage11.forcefield_config import validate_forcefield_stack


def test_default_forcefield_stack_has_no_warning():
    config = {"forcefield":{"ligand_parameterization_backend":"amber_gaff2_acpype","protein_ff":"amber19sb","ligand_ff":"gaff2","water_model":"opc"}}
    assert validate_forcefield_stack(config) == []


def test_gaff_with_charmm_stack_warns():
    config = {"forcefield":{"ligand_parameterization_backend":"amber_gaff2_acpype","protein_ff":"CHARMM36m","ligand_ff":"gaff2","water_model":"tip3p"}}
    assert validate_forcefield_stack(config)
