from __future__ import annotations

from syndesis.stage3.pose_labeling import label_one


def test_strict_native_like_requires_rmsd_and_sanity_pass() -> None:
    label = label_one(1.5, "pass", 2.0, 3.0)
    assert label["strict"] is True
    assert label["label"] == "strict_native_like"


def test_invalid_pose_not_native_like() -> None:
    label = label_one(1.0, "failed", 2.0, 3.0)
    assert label["invalid"] is True
    assert label["strict"] is False

