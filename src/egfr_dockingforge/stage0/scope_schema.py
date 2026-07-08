from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any

import yaml


REQUIRED_TOP_LEVEL_KEYS = (
    "project",
    "target",
    "ligand_scope",
    "receptor_scope",
    "validation_scope",
    "non_claims",
)

MAJOR_SCOPE_DECISION_KEYS = REQUIRED_TOP_LEVEL_KEYS


@dataclass(frozen=True)
class ValidationResult:
    valid: bool
    errors: list[str]


def load_yaml_mapping(path: str | Path) -> dict[str, Any]:
    with Path(path).open("r", encoding="utf-8") as handle:
        payload = yaml.safe_load(handle) or {}
    if not isinstance(payload, dict):
        raise ValueError(f"Expected YAML mapping in {path}")
    return payload


def validate_required_top_level_keys(scope: dict[str, Any]) -> list[str]:
    missing = [key for key in REQUIRED_TOP_LEVEL_KEYS if key not in scope]
    return [f"Missing required top-level key: {key}" for key in missing]


def validate_source_records(sources: dict[str, Any]) -> list[str]:
    errors: list[str] = []
    for source_id, record in sources.items():
        if not isinstance(record, dict):
            errors.append(f"Source `{source_id}` must be a mapping.")
            continue
        for field in ("supported_claim", "title"):
            if not record.get(field):
                errors.append(f"Source `{source_id}` missing required field `{field}`.")
        if not record.get("doi") and not record.get("url"):
            errors.append(f"Source `{source_id}` must include either `doi` or `url`.")
        if record.get("confidence") not in {"high", "medium", "low"}:
            errors.append(f"Source `{source_id}` confidence must be high, medium, or low.")
    return errors


def validate_decision_sources(scope: dict[str, Any], sources: dict[str, Any]) -> list[str]:
    known_sources = set(sources)
    errors: list[str] = []
    for key in MAJOR_SCOPE_DECISION_KEYS:
        section = scope.get(key)
        if not isinstance(section, dict):
            continue
        section_sources = section.get("decision_sources")
        if not isinstance(section_sources, list) or not section_sources:
            errors.append(f"Scope decision `{key}` must define a non-empty decision_sources list.")
            continue
        unknown = [source_id for source_id in section_sources if source_id not in known_sources]
        if unknown:
            errors.append(f"Scope decision `{key}` references unknown source ids: {unknown}")
    return errors


def validate_scope_policy(scope: dict[str, Any]) -> list[str]:
    errors: list[str] = []
    ligand_scope = scope.get("ligand_scope", {})
    receptor_scope = scope.get("receptor_scope", {})
    non_claims = scope.get("non_claims", {})

    excluded_ligands = " ".join(str(item).lower() for item in ligand_scope.get("excluded", []))
    if "covalent" not in excluded_ligands:
        errors.append("ligand_scope.excluded must exclude covalent inhibitors in v1.")

    mutation_policy = receptor_scope.get("mutation_policy", {})
    if mutation_policy.get("use_mutation_status_as") != "metadata_only":
        errors.append("receptor_scope.mutation_policy.use_mutation_status_as must be metadata_only.")
    if mutation_policy.get("selectivity_claim_allowed_v1") is not False:
        errors.append("receptor_scope.mutation_policy.selectivity_claim_allowed_v1 must be false.")

    if non_claims.get("no_experimental_activity_claim") is not True:
        errors.append("non_claims.no_experimental_activity_claim must be true.")
    if non_claims.get("no_experimentally_confirmed_inhibitor_discovery") is not True:
        errors.append("non_claims.no_experimentally_confirmed_inhibitor_discovery must be true.")

    return errors


def validate_scope(scope: dict[str, Any], sources: dict[str, Any]) -> ValidationResult:
    errors: list[str] = []
    errors.extend(validate_required_top_level_keys(scope))
    errors.extend(validate_source_records(sources))
    errors.extend(validate_decision_sources(scope, sources))
    errors.extend(validate_scope_policy(scope))
    return ValidationResult(valid=not errors, errors=errors)


def validate_scope_files(scope_path: str | Path, sources_path: str | Path) -> ValidationResult:
    scope = load_yaml_mapping(scope_path)
    sources = load_yaml_mapping(sources_path)
    return validate_scope(scope, sources)

