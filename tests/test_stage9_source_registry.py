from pathlib import Path

import yaml


def test_stage9_source_registry_has_required_entries_and_traceability():
    payload = yaml.safe_load(Path("data/references/stage9_sources.yaml").read_text())
    entries = payload["sources"]
    ids = {entry["source_id"] for entry in entries}
    required = {
        "toolmol_2026",
        "mollingo_2026",
        "mtmol_2025",
        "molclaw_2026",
        "reinvent4_2024",
        "reinvent_original_2017",
        "brics_rdkit",
        "matched_molecular_pairs",
        "matched_molecular_series",
        "qwen3_2025",
        "deepseek_r1_2025",
        "qwen3_coder_next_2026",
        "ollama_runtime",
        "rdkit",
        "pains_2010",
        "lipinski_2001",
        "veber_2002",
    }
    assert required.issubset(ids)
    assert all(entry.get("DOI") or entry.get("URL") for entry in entries)
