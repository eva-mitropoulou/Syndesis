from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

from egfr_dockingforge.common.io import write_table


def compute_score_hacking_cases(master: pd.DataFrame, paths: dict[str, Path]) -> pd.DataFrame:
    rows = []
    for row in master[master["score_hacking_flag"]].to_dict("records"):
        if row.get("delta_gnina_cnnscore", 0) > 0 and not row.get("binding_mode_preserved_flag", False):
            kind = "gnina_improves_binding_mode_breaks"
            worsened = "binding_mode_preserved_flag"
        elif row.get("delta_gnina_cnnscore", 0) > 0 and row.get("delta_pose_confidence", 0) < 0:
            kind = "gnina_improves_pose_confidence_worsens"
            worsened = "pose_confidence"
        elif row.get("delta_candidate_score", 0) > 0 and row.get("medchem_risk_score", 0) > 0.4:
            kind = "candidate_score_improves_medchem_worsens"
            worsened = "medchem_risk_score"
        else:
            kind = "score_improves_evidence_worsens"
            worsened = "structural_or_medchem_evidence"
        rows.append(
            {
                "analog_id": row["analog_id"],
                "strategy_id": row["strategy_id"],
                "seed_id": row["seed_id"],
                "score_hacking_type": kind,
                "improved_metric": "delta_gnina_cnnscore" if row.get("delta_gnina_cnnscore", 0) > 0 else "delta_candidate_score",
                "worsened_metric": worsened,
                "parent_value": 0.0,
                "analog_value": row.get("delta_gnina_cnnscore", row.get("delta_candidate_score", 0.0)),
                "severity": "high" if not row.get("binding_mode_preserved_flag", False) else "medium",
                "evidence_json": json.dumps({"delta_pose_confidence": row.get("delta_pose_confidence"), "delta_key_interaction_recall": row.get("delta_key_interaction_recall")}),
                "warnings_json": json.dumps([]),
            }
        )
    out = pd.DataFrame(rows, columns=["analog_id","strategy_id","seed_id","score_hacking_type","improved_metric","worsened_metric","parent_value","analog_value","severity","evidence_json","warnings_json"])
    write_table(paths["processed"] / "score_hacking_cases.parquet", out)
    write_table(paths["processed"] / "score_hacking_cases.csv", out)
    return out
