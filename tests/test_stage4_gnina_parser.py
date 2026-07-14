from __future__ import annotations

from pathlib import Path

import pytest

from syndesis.stage4.gnina_runner import _prepare_ligand_for_gnina
from syndesis.stage4.score_parser import parse_gnina_output
from syndesis.stage4.scoring_engines import parse_gnina_version


def test_gnina_parser_extracts_core_scores() -> None:
    text = """
    Affinity: -8.2
    CNNscore: 0.742
    CNNaffinity: 6.13
    """
    parsed = parse_gnina_output(text)
    assert parsed["gnina_empirical_affinity"] == -8.2
    assert parsed["cnnscore"] == 0.742
    assert parsed["cnnaffinity"] == 6.13
    assert parsed["cnn_vs"] is not None


def test_gnina_parser_handles_missing_cnn_vs() -> None:
    parsed = parse_gnina_output("CNNscore: 0.2\n")
    assert parsed["cnnscore"] == 0.2
    assert parsed["cnnaffinity"] is None
    assert parsed["cnn_vs"] is None


def test_gnina_version_parser() -> None:
    assert parse_gnina_version("gnina v1.3.3 master") == "gnina v1.3.3"


def test_gnina_pose_conversion_requires_openbabel(tmp_path: Path) -> None:
    pose = tmp_path / "pose.pdbqt"
    pose.write_text("END\n", encoding="utf-8")
    with pytest.raises(RuntimeError, match="Open Babel is required"):
        _prepare_ligand_for_gnina(pose, {"processed": tmp_path})
