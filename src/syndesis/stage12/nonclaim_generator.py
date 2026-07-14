from __future__ import annotations

REQUIRED_NON_CLAIMS = [
    "These candidates are computationally prioritized molecules, not experimentally confirmed EGFR inhibitors.",
    "Scores and poses are hypotheses requiring experimental validation.",
    "No claim is made about cellular potency, selectivity, pharmacokinetics, toxicity, or clinical relevance.",
    "Covalent EGFR inhibition is outside v1 scope.",
]

FORBIDDEN_PHRASES = [
    "confirmed inhibitor",
    "active compound",
    "validated hit",
    "clinically relevant",
    "selective EGFR inhibitor",
]


def nonclaim_text() -> str:
    return " ".join(REQUIRED_NON_CLAIMS)


def nonclaim_dict() -> dict[str, str]:
    return {
        "no_experimental_activity_claim": REQUIRED_NON_CLAIMS[0],
        "no_clinical_claim": REQUIRED_NON_CLAIMS[2],
        "no_selectivity_claim": REQUIRED_NON_CLAIMS[2],
        "no_covalent_claim": REQUIRED_NON_CLAIMS[3],
    }


def assert_nonclaims(text: str) -> None:
    lowered = text.lower()
    missing = [claim for claim in REQUIRED_NON_CLAIMS if claim not in text]
    forbidden = [phrase for phrase in FORBIDDEN_PHRASES if phrase in lowered]
    if missing or forbidden:
        raise ValueError(f"Non-claim policy failed; missing={missing}; forbidden={forbidden}")
