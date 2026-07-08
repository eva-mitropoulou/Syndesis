from __future__ import annotations

from pathlib import Path

import yaml


def test_stage6_source_registry_entries_have_doi_or_url() -> None:
    payload = yaml.safe_load(Path("data/references/stage6_sources.yaml").read_text(encoding="utf-8"))
    assert payload["sources"]
    for source in payload["sources"]:
        assert source.get("doi") or source.get("url")
