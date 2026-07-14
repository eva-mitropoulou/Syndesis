from __future__ import annotations

import json
import requests


def fetch_bindingdb_uniprot(uniprot_id: str, cutoff_nm: int, timeout_seconds: int, max_bytes: int = 5_000_000) -> list[dict]:
    url = f"https://bindingdb.org/rest/getLigandsByUniprots?uniprot={uniprot_id}&cutoff={cutoff_nm}&response=application/json"
    response = requests.get(url, timeout=timeout_seconds, stream=True)
    response.raise_for_status()
    chunks: list[bytes] = []
    total = 0
    for chunk in response.iter_content(chunk_size=65536):
        if not chunk:
            continue
        total += len(chunk)
        if total > max_bytes:
            raise RuntimeError(f"BindingDB response exceeded configured byte cap ({max_bytes}); use a pinned bulk download for full import.")
        chunks.append(chunk)
    text = b"".join(chunks).decode("utf-8").strip()
    return json.loads(text) if text else []
