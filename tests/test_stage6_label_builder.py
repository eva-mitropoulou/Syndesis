from __future__ import annotations

from syndesis.stage6.label_builder import RELEVANCE


def test_rank_relevance_values_are_valid() -> None:
    assert set(RELEVANCE.values()) <= {0, 1, 2, 3}
    assert RELEVANCE["high_confidence_native_like"] == 3
    assert RELEVANCE["wrong_binding_mode"] == 0
