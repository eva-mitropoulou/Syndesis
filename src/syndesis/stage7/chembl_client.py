from __future__ import annotations

import requests


BASE = "https://www.ebi.ac.uk/chembl/api/data"


def discover_egfr_target(uniprot_id: str) -> str:
    url = f"{BASE}/target?target_components__accession={uniprot_id}&target_type=SINGLE%20PROTEIN&format=json&limit=10"
    payload = requests.get(url, timeout=60).json()
    targets = payload.get("targets", [])
    if not targets:
        raise RuntimeError(f"ChEMBL returned no SINGLE PROTEIN target for {uniprot_id}")
    return targets[0]["target_chembl_id"]


def fetch_activities(target_chembl_id: str, limit: int) -> list[dict]:
    records: list[dict] = []
    offset = 0
    page = min(1000, limit)
    while len(records) < limit:
        url = f"{BASE}/activity?target_chembl_id={target_chembl_id}&standard_type__in=IC50,Ki,Kd,EC50&format=json&limit={page}&offset={offset}"
        payload = requests.get(url, timeout=90).json()
        batch = payload.get("activities", [])
        records.extend(batch)
        if not payload.get("page_meta", {}).get("next") or not batch:
            break
        offset += page
    return records[:limit]
