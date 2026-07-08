from __future__ import annotations

import math
from typing import Any


def to_nm(value: Any, unit: str | None) -> float | None:
    if value in (None, ""):
        return None
    try:
        val = float(value)
    except (TypeError, ValueError):
        return None
    unit_l = (unit or "nM").lower()
    if unit_l in {"nm", "nanomolar"}:
        return val
    if unit_l in {"um", "µm", "micromolar"}:
        return val * 1000.0
    if unit_l in {"m", "molar"}:
        return val * 1e9
    return None


def p_activity_from_nm(value_nm: float | None) -> float | None:
    if value_nm is None or value_nm <= 0:
        return None
    return -math.log10(value_nm * 1e-9)
