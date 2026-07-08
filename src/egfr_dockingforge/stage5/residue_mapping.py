from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pandas as pd

from egfr_dockingforge.common.io import write_table
from egfr_dockingforge.stage5.schemas import RESIDUE_MAP_COLUMNS

EGFR_RESIDUE_ROLES = {
    745: "catalytic_lys",
    762: "catalytic_glu",
    790: "gatekeeper",
    793: "hinge",
    797: "covalent_cys",
    855: "dfg_region",
    856: "dfg_region",
    857: "dfg_region",
}


def residue_role(uniprot_residue_number: Any) -> str:
    try:
        return EGFR_RESIDUE_ROLES.get(int(float(uniprot_residue_number)), "pocket_residue")
    except (TypeError, ValueError):
        return "unknown"


def normalize_interaction_key(row: pd.Series | dict[str, Any]) -> str:
    get = row.get if isinstance(row, dict) else row.get
    residue = get("uniprot_residue_number")
    interaction = get("interaction_type")
    role = get("residue_role") or residue_role(residue)
    if pd.isna(residue):
        residue = get("auth_seq_id")
    return f"{role}:{int(float(residue)) if pd.notna(residue) else 'NA'}:{interaction}"


def build_interaction_residue_map(inputs: dict[str, pd.DataFrame], paths: dict[str, Path]) -> pd.DataFrame:
    mapping = inputs["pocket_mapping"].copy()
    rows: list[dict[str, Any]] = []
    for _, row in mapping.iterrows():
        warnings: list[str] = []
        klifs = row.get("klifs_position")
        uniprot = row.get("uniprot_residue_number")
        if pd.notna(klifs):
            source = "klifs"
            confidence = "high"
        elif pd.notna(uniprot):
            source = "uniprot_fallback"
            confidence = "medium"
        else:
            source = "author_residue_fallback"
            confidence = "low"
            warnings.append("missing_uniprot_and_klifs_mapping")
        rows.append(
            {
                "receptor_id": row.get("receptor_id"),
                "pdb_id": row.get("pdb_id"),
                "auth_asym_id": row.get("auth_asym_id"),
                "residue_name": row.get("residue_name"),
                "auth_seq_id": row.get("auth_seq_id"),
                "uniprot_residue_number": uniprot,
                "klifs_position": klifs,
                "residue_role": residue_role(uniprot),
                "mapping_source": source,
                "mapping_confidence": confidence,
                "warnings_json": json.dumps(warnings),
            }
        )
    frame = pd.DataFrame(rows, columns=RESIDUE_MAP_COLUMNS)
    write_table(paths["processed"] / "interaction_residue_map.parquet", frame)
    write_table(paths["processed"] / "interaction_residue_map.csv", frame)
    return frame
