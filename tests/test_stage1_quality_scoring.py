from __future__ import annotations

from syndesis.stage1.quality_scoring import assign_quality
from syndesis.stage0.scope_schema import load_yaml_mapping


def test_tier_a_quality_assignment_for_good_complex() -> None:
    config = load_yaml_mapping("configs/stage1_cocrystal_benchmark.yaml")
    record = {
        "resolution_angstrom": 2.0,
        "ligand_heavy_atom_count": 25,
        "ligand_occupancy_mean": 1.0,
        "active_site_completeness_score": 1.0,
        "ligand_altloc_flag": False,
        "validation_report_available": True,
        "hard_exclusion_reasons": [],
    }
    quality = assign_quality(record, config)
    assert quality["quality_tier"] == "Tier A"
    assert quality["include_in_stage1_benchmark"] is True
    assert quality["quality_score"] is not None


def test_rejected_quality_has_reason_and_no_score() -> None:
    config = load_yaml_mapping("configs/stage1_cocrystal_benchmark.yaml")
    record = {
        "resolution_angstrom": 2.0,
        "ligand_heavy_atom_count": 25,
        "ligand_occupancy_mean": 1.0,
        "active_site_completeness_score": 1.0,
        "ligand_altloc_flag": False,
        "validation_report_available": True,
        "hard_exclusion_reasons": ["Covalent or likely covalent ligand."],
    }
    quality = assign_quality(record, config)
    assert quality["quality_tier"] == "Rejected"
    assert quality["include_in_stage1_benchmark"] is False
    assert quality["quality_score"] is None
    assert quality["exclusion_reason"]

