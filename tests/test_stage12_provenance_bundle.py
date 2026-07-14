from __future__ import annotations

import json
from pathlib import Path

import pytest

from syndesis.stage12.candidate_dossiers import run_stage12_all


CONFIG = Path("configs/stage12_candidate_dossiers.yaml")
pytestmark = pytest.mark.integration


def test_provenance_manifest_contains_required_sections() -> None:
    run_stage12_all(CONFIG)
    path = Path("reports/reproducibility_bundle/manifest.json")
    assert path.exists()
    manifest = json.loads(path.read_text(encoding="utf-8"))
    for key in ["project_name", "run_id", "git_commit", "environment", "config_hashes", "input_table_hashes", "output_table_hashes", "model_artifacts", "report_files", "tool_versions", "generated_at", "limitations"]:
        assert key in manifest
    assert manifest["config_hashes"]
    assert manifest["input_table_hashes"]
    assert manifest["output_table_hashes"]
