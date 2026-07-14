import json

import pandas as pd

from syndesis.stage10.ablation_tables import _stage11_md_ablation_row


def test_stage11_md_ablation_uses_completed_labels():
    row = _stage11_md_ablation_row(
        {
            "stage11_md_metrics": pd.DataFrame(
                {
                    "md_candidate_id": ["mdcand_001", "mdcand_002"],
                    "trajectory_analysis_status": ["complete", "complete"],
                }
            ),
            "stage11_md_pose_stability_labels": pd.DataFrame(
                {
                    "md_candidate_id": ["mdcand_001", "mdcand_002"],
                    "md_acceptance_flag": [False, False],
                }
            ),
        }
    )

    assert row["conclusion"] == "stage11_completed_no_candidates_passed_md_filter"
    # No stage9 acceptance table is provided, so the pre-MD rate cannot be derived
    # and must default to 0.0 (not the fabricated 1.0), yielding no fabricated effect.
    assert row["accepted_rate_change"] == 0.0
    assert json.loads(row["warnings_json"]) == ["pre_md_accept_rate_unavailable_defaulted_to_zero"]


def test_stage11_md_ablation_uses_actual_pre_md_rate():
    row = _stage11_md_ablation_row(
        {
            "stage11_md_metrics": pd.DataFrame(
                {
                    "md_candidate_id": ["mdcand_001", "mdcand_002"],
                    "trajectory_analysis_status": ["complete", "complete"],
                }
            ),
            "stage11_md_pose_stability_labels": pd.DataFrame(
                {
                    "md_candidate_id": ["mdcand_001", "mdcand_002"],
                    "md_acceptance_flag": [False, False],
                }
            ),
            "stage9_analog_acceptance": pd.DataFrame(
                {
                    "analog_id": ["a1", "a2", "a3", "a4"],
                    "accepted_flag": [False, False, False, False],
                }
            ),
        }
    )

    assert row["conclusion"] == "stage11_completed_no_candidates_passed_md_filter"
    # True pre-MD accepted rate is 0/4 = 0.0, post-MD is 0/2 = 0.0, so no fabricated effect.
    assert row["accepted_rate_change"] == 0.0
    assert json.loads(row["warnings_json"]) == []
