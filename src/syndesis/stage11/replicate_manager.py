from __future__ import annotations


def replicate_ids(quick: bool = True, count: int = 1) -> list[str]:
    return [f"rep{i:02d}" for i in range(1, count + 1)]
