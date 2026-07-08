from __future__ import annotations

from pathlib import Path

import yaml


def test_stage5_source_registry_entries_have_doi_or_url() -> None:
    sources = yaml.safe_load(Path("data/references/stage5_sources.yaml").read_text(encoding="utf-8"))
    required = {"prolif_2021", "plip_2015", "klifs_database", "klifs_2021_nar", "interaction_recovery_2024", "posebusters_2024", "oddt_2015_optional", "ifp_original_methods_optional", "kinase_ifp_binding_modes_optional"}
    assert required.issubset(sources)
    for source_id, source in sources.items():
        assert source["source_id"] == source_id
        assert source.get("doi") or source.get("url")
