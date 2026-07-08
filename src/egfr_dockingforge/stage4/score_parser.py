from __future__ import annotations

import re
from typing import Any

FLOAT_PATTERN = r"([-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?)"

LABEL_PATTERNS = {
    "cnnscore": [re.compile(rf"\bCNNscore\b\s*[:=]\s*{FLOAT_PATTERN}", re.IGNORECASE)],
    "cnnaffinity": [re.compile(rf"\bCNNaffinity\b\s*[:=]\s*{FLOAT_PATTERN}", re.IGNORECASE)],
    "gnina_empirical_affinity": [re.compile(rf"(?<!CNN)\bAffinity\b\s*[:=]\s*{FLOAT_PATTERN}", re.IGNORECASE)],
}


def to_float(value: Any) -> float | None:
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def parse_score_table_line(text: str, parsed: dict[str, float | None]) -> None:
    header: list[str] | None = None
    for raw in text.splitlines():
        line = raw.strip()
        if not line:
            continue
        if line.startswith("#"):
            candidate = line.lstrip("#").strip()
            if "CNNscore" in candidate and "CNNaffinity" in candidate:
                header = candidate.split()
            continue
        if header is None:
            continue
        values = line.split()
        if len(values) < len(header):
            continue
        mapping = dict(zip(header, values, strict=False))
        if parsed["cnnscore"] is None:
            parsed["cnnscore"] = to_float(mapping.get("CNNscore"))
        if parsed["cnnaffinity"] is None:
            parsed["cnnaffinity"] = to_float(mapping.get("CNNaffinity"))
        if parsed["gnina_empirical_affinity"] is None:
            parsed["gnina_empirical_affinity"] = to_float(mapping.get("Affinity"))


def parse_gnina_output(text: str) -> dict[str, float | None]:
    parsed: dict[str, float | None] = {
        "gnina_empirical_affinity": None,
        "cnnscore": None,
        "cnnaffinity": None,
        "cnn_vs": None,
    }
    for key, patterns in LABEL_PATTERNS.items():
        for pattern in patterns:
            matches = pattern.findall(text)
            if matches:
                parsed[key] = to_float(matches[-1])
                break
    parse_score_table_line(text, parsed)
    if parsed["cnnscore"] is not None and parsed["cnnaffinity"] is not None:
        parsed["cnn_vs"] = parsed["cnnscore"] * parsed["cnnaffinity"]
    return parsed

