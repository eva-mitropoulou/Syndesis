from __future__ import annotations

from syndesis.stage0.scope_schema import load_yaml_mapping


def test_stage1_source_registry_entries_have_doi_or_url() -> None:
    sources = load_yaml_mapping("data/references/stage1_sources.yaml")
    required = {
        "rcsb_search_api",
        "rcsb_data_api",
        "rcsb_mmcif_downloads",
        "rcsb_validation_reports",
        "klifs_database",
        "klifs_2021_nar",
        "kincore_database",
        "kincore_2019_pnas",
        "egfr_1m17_erlotinib",
        "egfr_1xkk_lapatinib",
        "egfr_4zau_azd9291_exclusion",
        "posebusters_2023",
        "interaction_recovery_2024",
        "leakproof_pdbbind_2023",
    }
    assert required.issubset(sources)
    for source_id, source in sources.items():
        assert source["source_id"] == source_id
        assert source.get("supported_claim")
        assert source.get("title")
        assert source.get("doi") or source.get("url")
        assert source.get("confidence") in {"high", "medium", "low"}

