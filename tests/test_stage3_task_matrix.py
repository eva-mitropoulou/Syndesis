from __future__ import annotations

from syndesis.stage3.task_matrix import task_type


def test_task_type_distinguishes_redocking_and_crossdocking() -> None:
    assert task_type("active_like", "active_like", "r1", "r1") == ("redocking_native_receptor", True)
    assert task_type("active_like", "active_like", "r1", "r2") == ("crossdocking_same_state", True)
    assert task_type("active_like", "inactive_like", "r1", "r2") == ("crossdocking_other_state", False)

