from __future__ import annotations

from egfr_dockingforge.stage0.scope_schema import load_yaml_mapping


def test_stage3_source_registry_entries_have_doi_or_url() -> None:
    sources = load_yaml_mapping("data/references/stage3_sources.yaml")
    required = {
        "autodock_vina_original_2010",
        "autodock_vina_1_2_2021",
        "autodock_gpu",
        "unidock_if_used",
        "diffdock_2022_optional",
        "unimol_docking_v2_optional",
        "posebusters_2023",
        "interaction_recovery_2024",
        "fair_comparison_docking_2024",
        "litpcba_docking_benchmark_2026",
        "docking_crossdocking_benchmark_general",
    }
    assert required.issubset(sources)
    for source_id, source in sources.items():
        assert source["source_id"] == source_id
        assert source.get("doi") or source.get("url")

