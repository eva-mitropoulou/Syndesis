from __future__ import annotations

import gzip
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import requests

from egfr_dockingforge.common.io import ensure_dir, write_json


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def normalize_pdb_id(value: str) -> str:
    cleaned = "".join(char for char in str(value).strip().upper() if char.isalnum())
    if len(cleaned) != 4:
        raise ValueError(f"Invalid PDB ID: {value!r}")
    return cleaned


def build_search_payload(config: dict[str, Any]) -> dict[str, Any]:
    query = config["rcsb"]["query"]
    target = config["target"]
    nodes: list[dict[str, Any]] = [
        {
            "type": "terminal",
            "service": "text",
            "parameters": {
                "attribute": "rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers.database_accession",
                "operator": "exact_match",
                "value": target.get("uniprot_accession", "P00533"),
            },
        },
        {
            "type": "terminal",
            "service": "text",
            "parameters": {
                "attribute": "rcsb_entry_info.nonpolymer_entity_count",
                "operator": "greater",
                "value": 0,
            },
        },
    ]
    methods = query.get("experimental_methods") or ["X-RAY DIFFRACTION"]
    nodes.append(
        {
            "type": "terminal",
            "service": "text",
            "parameters": {"attribute": "exptl.method", "operator": "in", "value": methods},
        }
    )
    if query.get("max_resolution_angstrom") is not None:
        nodes.append(
            {
                "type": "terminal",
                "service": "text",
                "parameters": {
                    "attribute": "rcsb_entry_info.resolution_combined",
                    "operator": "less_or_equal",
                    "value": float(query["max_resolution_angstrom"]),
                },
            }
        )
    return {
        "query": {"type": "group", "logical_operator": "and", "nodes": nodes},
        "return_type": "entry",
        "request_options": {
            "paginate": {"start": 0, "rows": int(query.get("max_structures", 80))},
            "sort": [{"sort_by": "rcsb_entry_info.resolution_combined", "direction": "asc"}],
            "scoring_strategy": "combined",
        },
    }


def search_pdb_ids(config: dict[str, Any]) -> list[str]:
    query = config["rcsb"]["query"]
    manual_ids = query.get("manual_pdb_ids") or []
    if query.get("use_manual_pdb_ids", False) and manual_ids:
        ids = [normalize_pdb_id(pdb_id) for pdb_id in manual_ids]
    else:
        response = requests.post(config["rcsb"]["search_url"], json=build_search_payload(config), timeout=90)
        response.raise_for_status()
        ids = [normalize_pdb_id(row["identifier"]) for row in response.json().get("result_set", [])]
    for control_id in config.get("controls", {}).get("required_reference_pdb_ids", []):
        ids.append(normalize_pdb_id(control_id))
    return list(dict.fromkeys(ids))[: int(query.get("max_structures", 80))]


def request_timeout(config: dict[str, Any], optional: bool = False) -> int:
    downloads = config.get("rcsb", {}).get("downloads", {})
    key = "optional_timeout_seconds" if optional else "timeout_seconds"
    return int(downloads.get(key, 8 if optional else 30))


def fetch_json(url: str, out_path: Path, force: bool = False, timeout: int = 30) -> dict[str, Any]:
    ensure_dir(out_path.parent)
    if out_path.exists() and out_path.stat().st_size > 0 and not force:
        return json.loads(out_path.read_text(encoding="utf-8"))
    response = requests.get(url, timeout=timeout)
    response.raise_for_status()
    payload = response.json()
    out_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return payload


def fetch_binary(url: str, out_path: Path, force: bool = False, timeout: int = 30) -> bool:
    ensure_dir(out_path.parent)
    if out_path.exists() and out_path.stat().st_size > 0 and not force:
        return True
    response = requests.get(url, timeout=timeout)
    if response.status_code == 404:
        return False
    response.raise_for_status()
    out_path.write_bytes(response.content)
    return True


def fetch_mmcif(pdb_id: str, config: dict[str, Any], out_dir: Path, force: bool = False) -> Path:
    pdb_id = normalize_pdb_id(pdb_id)
    url = config["rcsb"]["mmcif_url_template"].format(pdb_id=pdb_id)
    path = out_dir / f"{pdb_id.lower()}.cif"
    fetch_binary(url, path, force=force, timeout=request_timeout(config))
    return path


def fetch_entry_metadata(pdb_id: str, config: dict[str, Any], out_dir: Path, force: bool = False) -> dict[str, Any]:
    pdb_id = normalize_pdb_id(pdb_id)
    url = f"{config['rcsb']['data_api_url']}/entry/{pdb_id}"
    return fetch_json(url, out_dir / f"{pdb_id.lower()}_entry.json", force=force, timeout=request_timeout(config))


def fetch_chem_comp(comp_id: str, config: dict[str, Any], out_dir: Path, force: bool = False) -> Path | None:
    comp_id = str(comp_id).upper()
    url = config["rcsb"]["ccd_url_template"].format(comp_id=comp_id)
    path = out_dir / f"{comp_id}.cif"
    return path if fetch_binary(url, path, force=force, timeout=request_timeout(config, optional=True)) else None


def validation_middle(pdb_id: str) -> str:
    pdb_id_lower = pdb_id.lower()
    return pdb_id_lower[1:3]


def fetch_validation_files(
    pdb_id: str, config: dict[str, Any], out_dir: Path, force: bool = False
) -> tuple[Path | None, Path | None]:
    pdb_id = normalize_pdb_id(pdb_id)
    if not config.get("rcsb", {}).get("downloads", {}).get("fetch_validation_reports", True):
        return None, None
    pdb_id_lower = pdb_id.lower()
    values = {"pdb_id_lower": pdb_id_lower, "middle": validation_middle(pdb_id)}
    xml_url = config["rcsb"]["validation_xml_url_template"].format(**values)
    pdf_url = config["rcsb"]["validation_pdf_url_template"].format(**values)
    gz_path = out_dir / f"{pdb_id_lower}_validation.xml.gz"
    pdf_path = out_dir / f"{pdb_id_lower}_validation.pdf"
    xml_path: Path | None = None
    try:
        if fetch_binary(xml_url, gz_path, force=force, timeout=request_timeout(config, optional=True)):
            xml_path = out_dir / f"{pdb_id_lower}_validation.xml"
            if force or not xml_path.exists():
                xml_path.write_bytes(gzip.decompress(gz_path.read_bytes()))
    except Exception:
        xml_path = None
    try:
        pdf = pdf_path if fetch_binary(pdf_url, pdf_path, force=force, timeout=request_timeout(config, optional=True)) else None
    except Exception:
        pdf = None
    return xml_path, pdf


def write_download_manifest(records: list[dict[str, Any]], out_path: Path) -> None:
    write_json(out_path, {"downloaded_at": utc_now(), "records": records})


def file_record(path: Path, source_url: str, parser_version: str) -> dict[str, Any]:
    return {
        "path": str(path),
        "source_url": source_url,
        "download_timestamp": utc_now(),
        "checksum_sha256": sha256_file(path) if path.exists() else None,
        "parser_version": parser_version,
        "code_version": None,
    }
