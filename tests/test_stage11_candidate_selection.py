import pandas as pd

from egfr_dockingforge.stage11.candidate_selection import select_md_candidates


def test_candidate_selection_uses_top_ranked(tmp_path):
    ranked = pd.DataFrame([{"molecule_id":"mol_a","source":"chembl","standard_smiles":"CC","best_screening_pose_id":"pose1","best_target_receptor_id":"rec1","final_candidate_score":1.0}])
    config = {"selection":{"max_quick_md_candidates":1,"max_replicate_md_candidates":1}}
    out = select_md_candidates({"stage8_ranked_candidates": ranked}, config, {"processed": tmp_path})
    assert len(out) == 1
    assert out["selected_for_quick_md"].item()
    assert out["ligand_file"].item().endswith("data/processed/stage7/prepared_ligands/prep_mol_a.sdf")
