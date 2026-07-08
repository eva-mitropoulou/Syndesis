from __future__ import annotations

from pathlib import Path

import pandas as pd

from egfr_dockingforge.stage12.structure_figures import _svg_text


def write_md_stability_distribution(selection: pd.DataFrame, out_dir: Path) -> Path:
    counts = selection["md_stability_label_if_available"].fillna("not_available").value_counts().to_dict()
    return _svg_text(out_dir / "md_stability_distribution.svg", "MD stability distribution", [f"{key}: {value}" for key, value in counts.items()])
