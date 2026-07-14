from syndesis.stage11.mdp_templates import production_mdp


def test_production_mdp_is_npt_and_unrestrained_ligand():
    config = {"md":{"replicate_production_ns":100,"quick_production_ns":20,"timestep_fs":2,"temperature_k":300,"pressure_bar":1.0}}
    text = production_mdp(config)
    assert "pcoupl = Parrinello-Rahman" in text
    assert "no ligand restraints" in text
    assert "define = -DPOSRES" not in text
