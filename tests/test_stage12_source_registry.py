from __future__ import annotations

from pathlib import Path

import yaml


REQUIRED = {
    "fair_principles_2016", "model_cards_2019", "datasheets_for_datasets_2018",
    "joss_research_software_review", "jupyter_reproducibility_biomedical_2023",
    "prolif_2021", "plip_2015", "molstar_2021", "gnina_2021", "posebusters_2024",
    "interaction_recovery_2024", "ultralarge_docking_2019", "rdkit", "shap_2017",
}


def test_stage12_source_registry_has_required_entries_and_links() -> None:
    payload = yaml.safe_load(Path("data/references/stage12_sources.yaml").read_text(encoding="utf-8"))
    sources = payload["sources"]
    by_id = {source["source_id"]: source for source in sources}
    assert REQUIRED.issubset(by_id)
    for source in sources:
        assert source.get("DOI") or source.get("URL")
        assert source.get("used_for")
        assert source.get("confidence") in {"high", "medium", "low"}
