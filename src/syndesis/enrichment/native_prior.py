from __future__ import annotations

import json
from collections.abc import Iterable

import pandas as pd


def parse_fingerprint(value: str | Iterable[str]) -> set[str]:
    if isinstance(value, str):
        value = json.loads(value) if value else []
    return {str(bit) for bit in value}


def native_union(
    native_fingerprints: pd.DataFrame,
    excluded_receptors: Iterable[str] = (),
) -> tuple[set[str], list[str]]:
    excluded = {str(value).lower() for value in excluded_receptors}
    selected = native_fingerprints[
        ~native_fingerprints["receptor_id"].astype(str).str.lower().isin(excluded)
    ]
    if selected.empty:
        raise ValueError("The native interaction prior contains no included complexes.")
    union = set().union(
        *(parse_fingerprint(value) for value in selected["fingerprint_sparse_json"])
    )
    if not union:
        raise ValueError("The native interaction prior contains no interaction bits.")
    return union, selected["receptor_id"].astype(str).tolist()


def recall(observed: set[str], target: set[str]) -> float:
    if not target:
        raise ValueError("Interaction recall requires a non-empty target prior.")
    return len(observed & target) / len(target)


def jaccard(observed: set[str], target: set[str]) -> float:
    union = observed | target
    return len(observed & target) / len(union) if union else 0.0
