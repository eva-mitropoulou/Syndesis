from __future__ import annotations

from syndesis.stage0.scope_schema import load_yaml_mapping


def test_stage2_source_registry_entries_have_doi_or_url() -> None:
    sources = load_yaml_mapping("data/references/stage2_sources.yaml")
    required = {
        "klifs_database",
        "klifs_2021_nar",
        "kincore_database",
        "kincore_2019_pnas",
        "kincore_alignment_2019_scirep",
        "egfr_1m17_erlotinib_active_like",
        "egfr_1xkk_lapatinib_inactive_like",
        "receptor_flexibility_ensemble_docking_review",
        "receptor_conformational_ensembles_virtual_screening",
        "posebusters_2023",
        "interaction_recovery_2024",
    }
    assert required.issubset(sources)
    for source_id, source in sources.items():
        assert source["source_id"] == source_id
        assert source.get("supported_claim")
        assert source.get("title")
        assert source.get("doi") or source.get("url")
        assert source.get("confidence") in {"high", "medium", "low"}

