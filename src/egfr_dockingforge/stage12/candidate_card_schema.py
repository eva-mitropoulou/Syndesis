from __future__ import annotations

REQUIRED_CARD_FIELDS = [
    "card_version", "final_candidate_id", "molecule_id", "source", "subsource", "screening_role",
    "standard_smiles", "inchi_key", "scaffold_id", "novelty_bucket", "closest_known_egfr_ligand",
    "parent_analog_lineage", "best_pose", "scores", "interactions", "md", "medchem",
    "evidence_summary", "non_claims", "provenance",
]


def validate_candidate_card(card: dict) -> None:
    missing = [field for field in REQUIRED_CARD_FIELDS if field not in card]
    if missing:
        raise ValueError(f"Missing candidate card fields: {missing}")
    for nested in ["closest_known_egfr_ligand", "parent_analog_lineage", "best_pose", "scores", "interactions", "md", "medchem", "evidence_summary", "non_claims", "provenance"]:
        if not isinstance(card[nested], dict):
            raise TypeError(f"Candidate card field {nested} must be an object.")
