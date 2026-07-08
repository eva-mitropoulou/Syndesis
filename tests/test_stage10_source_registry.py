from pathlib import Path

import yaml


def test_stage10_source_registry_entries_trace_claims():
    payload = yaml.safe_load(Path("data/references/stage10_sources.yaml").read_text())
    entries = payload["sources"]
    required = {
        "guacamol_2019",
        "moses_2018",
        "frechet_chemnet_distance_2018",
        "dockstring_2021",
        "sbdd_benchmark_2024",
        "genbench3d_2024",
        "litpcba_rescoring_supervised_ml_2026",
        "posebusters_2024",
        "interaction_recovery_2024",
        "reinvent4_2024",
        "toolmol_2026",
        "mollingo_2026",
        "scipy_statistics_docs",
        "bootstrap_confidence_intervals_reference",
        "benjamini_hochberg_1995",
    }
    assert required.issubset({entry["source_id"] for entry in entries})
    assert all(entry.get("DOI") or entry.get("URL") for entry in entries)
