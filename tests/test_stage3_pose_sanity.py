from __future__ import annotations

import pandas as pd

from syndesis.stage3.pose_sanity import sanity_for_poses


def test_missing_pose_file_triggers_atom_loss(tmp_path) -> None:
    poses = pd.DataFrame([{"pose_id": "p1", "pose_file": str(tmp_path / "missing.pdb")}])
    sanity = sanity_for_poses(poses, {"processed": tmp_path})
    assert sanity.loc[0, "atom_loss_flag"] is True or sanity.loc[0, "atom_loss_flag"] == True
    assert sanity.loc[0, "sanity_status"] == "failed"

