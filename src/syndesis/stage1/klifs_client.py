from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import requests

from syndesis.common.io import ensure_dir


NULL_KLIFS_METADATA = {
    "klifs_structure_id": None,
    "klifs_kinase_id": None,
    "klifs_ligand_id": None,
    "klifs_pocket_id": None,
    "klifs_dfg_state": None,
    "klifs_ac_helix_state": None,
}


def fetch_klifs_metadata(pdb_id: str, config: dict[str, Any], out_dir: Path) -> tuple[dict[str, Any], list[str]]:
    if not config.get("klifs", {}).get("enabled", True):
        return dict(NULL_KLIFS_METADATA), ["KLIFS metadata disabled."]
    ensure_dir(out_dir)
    pdb_id = pdb_id.upper()
    out_path = out_dir / f"{pdb_id.lower()}_klifs.json"
    warnings: list[str] = []
    if out_path.exists() and out_path.stat().st_size > 0:
        try:
            payload = json.loads(out_path.read_text(encoding="utf-8"))
            return normalize_klifs_payload(payload), warnings
        except Exception as exc:
            warnings.append(f"Could not parse cached KLIFS metadata for {pdb_id}: {exc}")

    # KLIFS deployments have changed endpoint names over time. Try conservative read-only
    # endpoints and fail softly rather than inventing kinase annotations.
    base_url = config.get("klifs", {}).get("base_url", "https://klifs.net/api").rstrip("/")
    timeout = int(config.get("klifs", {}).get("timeout_seconds", 5))
    candidate_urls = [
        f"{base_url}/structures_pdb_list?pdb-codes={pdb_id}",
        f"{base_url}/structures_pdb_list?pdb-codes={pdb_id.lower()}",
    ]
    for url in candidate_urls:
        try:
            response = requests.get(url, timeout=timeout)
            if response.status_code == 404:
                continue
            response.raise_for_status()
            payload = response.json()
            out_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
            return normalize_klifs_payload(payload), warnings
        except Exception as exc:
            warnings.append(f"KLIFS lookup failed at {url}: {exc}")
    return dict(NULL_KLIFS_METADATA), warnings or [f"KLIFS metadata unavailable for {pdb_id}."]


def normalize_klifs_payload(payload: Any) -> dict[str, Any]:
    row: dict[str, Any] = {}
    if isinstance(payload, list) and payload:
        first = payload[0]
        row = first if isinstance(first, dict) else {}
    elif isinstance(payload, dict):
        if isinstance(payload.get("data"), list) and payload["data"]:
            row = payload["data"][0]
        else:
            row = payload
    return {
        "klifs_structure_id": row.get("structure_ID") or row.get("structure_id") or row.get("id"),
        "klifs_kinase_id": row.get("kinase_ID") or row.get("kinase_id"),
        "klifs_ligand_id": row.get("ligand_ID") or row.get("ligand_id"),
        "klifs_pocket_id": row.get("pocket_ID") or row.get("pocket_id"),
        "klifs_dfg_state": row.get("DFG"),
        "klifs_ac_helix_state": row.get("aC_helix"),
    }
