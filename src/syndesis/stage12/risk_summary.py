from __future__ import annotations

import json
from typing import Any


def parse_json_list(value: Any) -> list:
    if value is None:
        return []
    if isinstance(value, list):
        return value
    if isinstance(value, float) and value != value:
        return []
    if isinstance(value, str):
        if not value:
            return []
        try:
            payload = json.loads(value)
            return payload if isinstance(payload, list) else [payload]
        except json.JSONDecodeError:
            return [value]
    return [value]


def merge_warnings(*values: Any) -> str:
    merged: list[str] = []
    for value in values:
        for item in parse_json_list(value):
            text = str(item)
            if text and text not in merged:
                merged.append(text)
    return json.dumps(merged)


def main_risks(row: dict) -> list[str]:
    risks = parse_json_list(row.get("risk_flags_json"))
    if row.get("decision_label") == "low_confidence_rejected":
        risks.append("low pose-confidence or weak interaction recovery")
    md_label = str(row.get("md_stability_label_if_available", "")).lower()
    md_decision = str(row.get("final_md_decision", "")).lower()
    if md_label.startswith("md_failed"):
        risks.append("no production MD stability evidence")
    if md_label == "md_unstable" or row.get("decision_label") == "md_unstable_rejected" or md_decision == "fail_md_stability":
        risks.append("MD pose instability (did not meet stability thresholds)")
    return sorted(set(str(risk) for risk in risks if risk))
