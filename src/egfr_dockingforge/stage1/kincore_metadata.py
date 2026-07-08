from __future__ import annotations

from typing import Any


KNOWN_STATE_HINTS: dict[str, dict[str, Any]] = {
    "1M17": {
        "kincore_state": "active-like",
        "dfg_state": "DFGin",
        "chelix_state": None,
        "saltbridge_state": None,
        "hrd_state": None,
        "activation_loop_state": "active-like",
    },
    "1XKK": {
        "kincore_state": "inactive-like",
        "dfg_state": None,
        "chelix_state": None,
        "saltbridge_state": None,
        "hrd_state": None,
        "activation_loop_state": "inactive-like",
    },
}


NULL_KINASE_STATE = {
    "kincore_state": None,
    "dfg_state": None,
    "chelix_state": None,
    "saltbridge_state": None,
    "hrd_state": None,
    "activation_loop_state": None,
}


def _clean_state(value: Any) -> str | None:
    if value in {None, "", "?", "."}:
        return None
    return str(value).strip().lower()


def infer_from_klifs(klifs_metadata: dict[str, Any] | None) -> tuple[dict[str, Any] | None, list[str]]:
    if not klifs_metadata:
        return None, []
    dfg = _clean_state(klifs_metadata.get("klifs_dfg_state"))
    ac_helix = _clean_state(klifs_metadata.get("klifs_ac_helix_state"))
    if not dfg and not ac_helix:
        return None, []

    dfg_state = "DFGout" if dfg == "out" else "DFGin" if dfg == "in" else None
    chelix_state = "C-helix-out" if ac_helix == "out" else "C-helix-in" if ac_helix == "in" else None
    if dfg_state == "DFGout":
        activity = "inactive-like"
        activation = "inactive-like"
    elif dfg_state == "DFGin" and chelix_state == "C-helix-in":
        activity = "active-like"
        activation = "active-like"
    elif dfg_state == "DFGin" and chelix_state == "C-helix-out":
        activity = "inactive-like"
        activation = "inactive-like"
    else:
        activity = None
        activation = None

    metadata = {
        **NULL_KINASE_STATE,
        "kincore_state": activity,
        "dfg_state": dfg_state,
        "chelix_state": chelix_state,
        "activation_loop_state": activation,
    }
    warnings = []
    if activity is None:
        warnings.append("KLIFS DFG/aC-helix state present but insufficient for active/inactive assignment.")
    return metadata, warnings


def infer_kincore_metadata(
    pdb_id: str,
    _chain_id: str,
    _config: dict[str, Any],
    klifs_metadata: dict[str, Any] | None = None,
) -> tuple[dict[str, Any], list[str]]:
    pdb_id = pdb_id.upper()
    if pdb_id in KNOWN_STATE_HINTS:
        return dict(KNOWN_STATE_HINTS[pdb_id]), []
    inferred, warnings = infer_from_klifs(klifs_metadata)
    if inferred is not None:
        return inferred, warnings
    return dict(NULL_KINASE_STATE), [f"KinCore/KLIFS state metadata unavailable for {pdb_id}; fields left null."]
