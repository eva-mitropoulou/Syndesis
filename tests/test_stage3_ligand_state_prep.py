from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest


PROJECT_ROOT = Path(__file__).resolve().parents[1]


def test_native_reference_files_are_separate_when_outputs_exist() -> None:
    path = PROJECT_ROOT / "data/processed/stage3/ligand_docking_prep.parquet"
    if not path.exists():
        pytest.skip("Stage 3 prep has not been generated.")
    frame = pd.read_parquet(path)
    assert (frame["native_ligand_file"] != frame["prepared_ligand_file"]).all()
    assert frame["immutable_reference_pose_file"].map(lambda p: Path(p).exists()).all()

