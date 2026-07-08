from pathlib import Path

import yaml


def test_stage11_source_registry_traceability():
    entries = yaml.safe_load(Path("data/references/stage11_sources.yaml").read_text())["sources"]
    required = {"charmm36m_2017","cgenff_2010","cgenff_program_2016","gromacs_reference","gromacs_2024_or_installed_version","charmm_tip3p","prolif_2021","protein_ligand_md_best_practices","md_reproducibility_reference","posebusters_2024","interaction_recovery_2024"}
    assert required.issubset({e["source_id"] for e in entries})
    assert all(e.get("DOI") or e.get("URL") for e in entries)
