from __future__ import annotations

import json
import pandas as pd


def aggregate_known(measurements: pd.DataFrame) -> pd.DataFrame:
    if measurements.empty:
        return pd.DataFrame()
    rows = []
    for (molecule_id, endpoint), group in measurements.groupby(["molecule_id", "endpoint_type"], dropna=False):
        vals = group["standard_value_nM"].dropna()
        rows.append(
            {
                "molecule_id": molecule_id,
                "standard_smiles": group["molecule_smiles_standard"].iloc[0],
                "inchi_key": group["inchi_key"].iloc[0],
                "source_list": ",".join(sorted(set(group["source"]))),
                "source_record_ids_json": json.dumps(group["source_record_id"].astype(str).tolist()),
                "known_activity_status": "known_measured",
                "endpoint_type": endpoint,
                "assay_context_group": "biochemical_or_binding",
                "median_p_activity": group["p_activity"].median(),
                "median_standard_value_nM": vals.median() if not vals.empty else None,
                "min_standard_value_nM": vals.min() if not vals.empty else None,
                "max_standard_value_nM": vals.max() if not vals.empty else None,
                "num_measurements": int(len(group)),
                "num_assays": int(group["assay_id"].nunique()),
                "num_documents": int(group["document_id"].nunique()),
                "conflict_flag": bool(len(vals) > 1 and vals.max() / max(vals.min(), 1e-9) > 100),
                "activity_confidence": "medium",
                "mutation_context_summary": "",
                "covalent_flag": False,
                "warhead_flag": False,
                "approved_or_clinical_flag": False,
                "native_ligand_flag": False,
                "include_as_known_reference": True,
                "include_as_positive_control": True,
                "include_in_retrospective_series": True,
                "exclusion_reason": "",
                "warnings_json": json.dumps([]),
            }
        )
    return pd.DataFrame(rows)
