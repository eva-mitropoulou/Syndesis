import pandas as pd

from egfr_dockingforge.stage11.stability_scoring import score_md_stability


def test_missing_md_labels_failed_setup(tmp_path):
    candidates = pd.DataFrame([{"md_candidate_id":"c1","molecule_id":"mol","source":"s","best_pose_id":"p","parent_molecule_id":"mol","stage10_strategy_name":""}])
    metrics = pd.DataFrame([{"md_candidate_id":"c1","md_system_id":"s1","replicate_id":"rep01"}])
    labels, summary, post = score_md_stability(candidates, metrics, pd.DataFrame(), {"processed": tmp_path})
    assert labels["md_stability_label"].item() == "md_failed_setup"
    assert not bool(post["accepted_post_md_flag"].item())
