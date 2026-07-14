from __future__ import annotations

import pandas as pd


def write_disabled_active_learning(paths: dict) -> pd.DataFrame:
    df = pd.DataFrame(columns=["iteration_id", "candidate_pool_size", "num_docked_this_iteration", "surrogate_model_type", "acquisition_function", "selected_molecule_ids_json", "top_score_recovery_estimate", "stopping_reason", "warnings_json"])
    df.to_parquet(paths["processed"] / "active_learning_iterations.parquet", index=False)
    return df
