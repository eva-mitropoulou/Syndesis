from __future__ import annotations

from pathlib import Path

import pandas as pd

from egfr_dockingforge.stage12.final_candidate_selection import build_final_candidate_table


CONFIG = Path("configs/stage12_candidate_dossiers.yaml")


def test_known_controls_are_separate_from_generated_candidates() -> None:
    build_final_candidate_table(CONFIG)
    selection = pd.read_parquet("data/processed/stage12/final_candidate_selection.parquet")
    known = selection[selection["screening_role"].eq("known_activity_reference")]
    generated = selection[selection["screening_role"].eq("generated_analog_negative_control")]
    assert not known.empty
    assert not generated.empty
    assert set(known["source"]) == {"chembl_known_ligand"}
    assert generated["source"].str.contains("stage9").all()


def test_md_missing_candidates_do_not_receive_md_stable_claims() -> None:
    build_final_candidate_table(CONFIG)
    selection = pd.read_parquet("data/processed/stage12/final_candidate_selection.parquet")
    label_col = selection["md_stability_label_if_available"].astype(str).str.lower()
    stable_claims = {"md_stable", "stable", "moderately_stable"}

    # A candidate may carry an md_stable claim ONLY if it actually has MD
    # evidence (non-null interaction persistence from a completed trajectory).
    # Candidates without MD evidence must be labelled "not_available".
    has_md_evidence = selection["md_key_interaction_persistence_if_available"].notna()
    stable_rows = label_col.isin(stable_claims)
    assert (stable_rows <= has_md_evidence).all(), "md_stable claim without MD evidence"
    missing = label_col.isin({"not_available", "nan", ""})
    assert not (missing & stable_rows).any(), "MD-missing candidate received a stability claim"

    # MD-unstable candidates must be rejected as such.
    unstable = label_col.eq("md_unstable")
    assert selection.loc[unstable, "decision_label"].eq("md_unstable_rejected").all()
    # MD-stable known controls must be recognised (not rejected as unstable).
    stable = label_col.eq("md_stable")
    assert not selection.loc[stable, "decision_label"].eq("md_unstable_rejected").any()
    assert set(selection["novelty_bucket"].dropna())
