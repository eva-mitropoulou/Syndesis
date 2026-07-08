from __future__ import annotations

from pathlib import Path

from egfr_dockingforge.stage0.scope_schema import (
    REQUIRED_TOP_LEVEL_KEYS,
    load_yaml_mapping,
    validate_scope_files,
)


PROJECT_ROOT = Path(__file__).resolve().parents[1]


def test_project_scope_required_yaml_keys_exist() -> None:
    scope = load_yaml_mapping(PROJECT_ROOT / "configs/project_scope.yaml")
    assert set(REQUIRED_TOP_LEVEL_KEYS).issubset(scope)


def test_project_scope_validates_against_sources() -> None:
    result = validate_scope_files(
        PROJECT_ROOT / "configs/project_scope.yaml",
        PROJECT_ROOT / "data/references/stage0_sources.yaml",
    )
    assert result.valid, result.errors


def test_non_claims_include_no_experimental_activity_claim() -> None:
    scope = load_yaml_mapping(PROJECT_ROOT / "configs/project_scope.yaml")
    assert scope["non_claims"]["no_experimental_activity_claim"] is True
    assert scope["non_claims"]["no_experimentally_confirmed_inhibitor_discovery"] is True


def test_covalent_inhibitors_are_excluded_in_v1() -> None:
    scope = load_yaml_mapping(PROJECT_ROOT / "configs/project_scope.yaml")
    excluded = " ".join(scope["ligand_scope"]["excluded"]).lower()
    assert "covalent inhibitors" in excluded
    assert "reversible covalent inhibitors" in excluded
    assert "irreversible covalent inhibitors" in excluded


def test_mutation_policy_is_metadata_only() -> None:
    scope = load_yaml_mapping(PROJECT_ROOT / "configs/project_scope.yaml")
    policy = scope["receptor_scope"]["mutation_policy"]
    assert policy["use_mutation_status_as"] == "metadata_only"
    assert policy["selectivity_claim_allowed_v1"] is False

