from __future__ import annotations

from syndesis.stage6.model_selection import _group_sizes


def test_ranker_group_sizes_preserve_order() -> None:
    assert _group_sizes(__import__("pandas").Series(["a", "a", "b", "b", "b"])) == [2, 3]
