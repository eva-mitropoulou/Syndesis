from __future__ import annotations

from typing import Any


def assign_quality(record: dict[str, Any], config: dict[str, Any]) -> dict[str, Any]:
    filters = config["filters"]
    reasons = list(record.get("hard_exclusion_reasons") or [])
    resolution = record.get("resolution_angstrom")
    heavy = record.get("ligand_heavy_atom_count") or 0
    occupancy = record.get("ligand_occupancy_mean")
    completeness = record.get("active_site_completeness_score") or 0.0

    if resolution is not None and resolution > float(filters.get("max_primary_resolution_angstrom", 3.5)):
        reasons.append("Resolution exceeds Stage 1 maximum.")
    if heavy < int(filters.get("min_ligand_heavy_atoms", 12)):
        reasons.append("Ligand heavy atom count below Stage 1 minimum.")
    if occupancy is not None and occupancy < float(filters.get("min_ligand_occupancy", 0.5)):
        reasons.append("Ligand occupancy below Stage 1 minimum.")

    include = not reasons
    score = 100.0
    if resolution is not None:
        score -= max(0.0, resolution - 1.0) * 8.0
    if occupancy is not None:
        score -= max(0.0, 1.0 - occupancy) * 20.0
    score -= max(0.0, 1.0 - completeness) * 25.0
    if record.get("ligand_altloc_flag"):
        score -= 10.0
    if not record.get("validation_report_available"):
        score -= 5.0
    score = round(max(0.0, min(100.0, score)), 2)

    if not include:
        tier = "Rejected"
        score = None
    elif resolution is not None and resolution <= float(filters.get("max_tier_a_resolution_angstrom", 2.5)) and completeness >= 0.875 and heavy >= 12:
        tier = "Tier A"
    elif resolution is not None and resolution <= float(filters.get("max_tier_b_resolution_angstrom", 3.0)) and completeness >= 0.75:
        tier = "Tier B"
    else:
        tier = "Tier C"

    return {
        "quality_tier": tier,
        "quality_score": score,
        "include_in_stage1_benchmark": include,
        "exclusion_reason": "; ".join(reasons) if reasons else None,
    }

