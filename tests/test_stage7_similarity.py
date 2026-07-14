from __future__ import annotations

from syndesis.stage7.similarity import novelty_bucket


def test_novelty_thresholds() -> None:
    assert novelty_bucket(1.0, True, True) == "known_duplicate"
    assert novelty_bucket(0.8, False, True) == "close_analog"
    assert novelty_bucket(0.5, False, False) == "medium_similarity"
    assert novelty_bucket(0.1, False, False) == "scaffold_novel"
