import pandas as pd

from syndesis.stage9.edit_site_detection import detect_edit_sites


def test_edit_sites_exist_and_protected_sites_have_no_allowed_classes(tmp_path):
    seeds = pd.DataFrame(
        [
            {
                "seed_id": "seed_001",
                "molecule_id": "mol1",
                "standard_smiles": "Cc1ccccc1",
            }
        ]
    )
    config = {"edit_sites": {"max_sites_per_seed": 6}, "transforms": {"enabled_classes": ["halogen_scan"]}}
    out = detect_edit_sites(seeds, config, {"processed": tmp_path})
    assert not out.empty
    protected = out[out["protected_region_flag"]]
    if not protected.empty:
        assert (protected["allowed_transformation_classes_json"] == "[]").all()
