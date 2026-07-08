from __future__ import annotations

from pathlib import Path

import yaml


def test_stage4_sources_have_doi_or_url() -> None:
    path = Path("data/references/stage4_sources.yaml")
    sources = yaml.safe_load(path.read_text(encoding="utf-8"))
    required = {"supported_claim", "title", "source_type", "used_for", "confidence"}
    for source_id, source in sources.items():
        assert source.get("source_id") == source_id
        assert required.issubset(source)
        assert source.get("doi") or source.get("url")

