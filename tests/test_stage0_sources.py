from __future__ import annotations

from pathlib import Path

from egfr_dockingforge.stage0.scope_schema import load_yaml_mapping


PROJECT_ROOT = Path(__file__).resolve().parents[1]


def test_every_stage0_source_has_required_metadata() -> None:
    sources = load_yaml_mapping(PROJECT_ROOT / "data/references/stage0_sources.yaml")
    assert sources
    for source_id, source in sources.items():
        assert source.get("supported_claim"), source_id
        assert source.get("title"), source_id
        assert source.get("doi") or source.get("url"), source_id


def test_scope_decision_source_ids_exist_in_registry() -> None:
    scope = load_yaml_mapping(PROJECT_ROOT / "configs/project_scope.yaml")
    sources = load_yaml_mapping(PROJECT_ROOT / "data/references/stage0_sources.yaml")
    known_source_ids = set(sources)
    for section_name in (
        "project",
        "target",
        "ligand_scope",
        "receptor_scope",
        "validation_scope",
        "non_claims",
    ):
        assert set(scope[section_name]["decision_sources"]).issubset(known_source_ids)

