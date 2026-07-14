from __future__ import annotations

import pandas as pd

from syndesis.stage8.source_aware_selection import rank_candidates


def test_source_aware_ranking_outputs_known_control_bucket(tmp_path) -> None:
    agg = pd.DataFrame({"molecule_id": ["m1"], "source": ["chembl_known_ligand"], "screening_role": ["known_activity_reference"], "novelty_bucket": ["close_analog"], "standard_smiles": ["CC"], "best_screening_pose_id": ["p1"], "best_target_receptor_id": ["r"], "best_receptor_state": ["active"], "final_candidate_score": [1.0], "best_pose_confidence": [0.8], "best_gnina_cnnscore": [0.7], "best_key_interaction_recall_consensus": [0.5], "best_ifp_tanimoto_to_consensus": [0.4], "medchem_risk_score": [0], "tanimoto_to_closest_known": [1.0], "closest_known_molecule_id": ["m1"], "scaffold_id": ["s"], "candidate_decision_label": ["accepted_candidate_control"], "candidate_decision_reason": ["accepted"], "medchem_flags_json": ["[]"], "warnings_json": ["[]"]})
    ranked, _ = rank_candidates(agg, {"processed": tmp_path})
    assert ranked["candidate_bucket"].iloc[0] == "known_control_recovered"
