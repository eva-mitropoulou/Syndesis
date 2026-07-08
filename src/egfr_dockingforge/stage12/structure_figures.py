from __future__ import annotations

from pathlib import Path

import pandas as pd


def _svg_text(path: Path, title: str, lines: list[str], width: int = 720, height: int = 420) -> Path:
    body = [f'<text x="24" y="42" font-size="22" font-family="Arial" font-weight="700">{title}</text>']
    y = 82
    for line in lines:
        escaped = str(line).replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")
        body.append(f'<text x="24" y="{y}" font-size="16" font-family="Arial">{escaped}</text>')
        y += 28
    path.write_text(
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">'
        '<rect width="100%" height="100%" fill="#f8fafc"/>'
        '<rect x="12" y="12" width="696" height="396" fill="#ffffff" stroke="#94a3b8"/>'
        + "".join(body)
        + "</svg>\n",
        encoding="utf-8",
    )
    return path


def write_candidate_structure_svg(row: dict, out_dir: Path) -> Path:
    try:
        from rdkit import Chem
        from rdkit.Chem import Draw

        mol = Chem.MolFromSmiles(row["standard_smiles"])
        if mol is not None:
            path = out_dir / f"{row['final_candidate_id']}_2d.svg"
            path.write_text(Draw.MolsToGridImage([mol], molsPerRow=1, subImgSize=(380, 280), useSVG=True), encoding="utf-8")
            return path
    except Exception:
        pass
    return _svg_text(out_dir / f"{row['final_candidate_id']}_2d.svg", row["final_candidate_id"], [row.get("standard_smiles", "structure unavailable")])


def write_top_candidate_grid(selection: pd.DataFrame, out_dir: Path) -> Path:
    return _svg_text(
        out_dir / "top_candidate_2d_grid.svg",
        "Top candidate 2D grid",
        [f"{r.final_candidate_id}: {r.molecule_id}" for r in selection.head(12).itertuples()],
        width=960,
        height=540,
    )


def write_pose_panels(selection: pd.DataFrame, out_dir: Path) -> Path:
    return _svg_text(
        out_dir / "top_candidate_3d_pose_panels.svg",
        "Top candidate pose panels",
        [f"{r.final_candidate_id}: {r.best_receptor_state or 'pose missing'} / {r.best_pose_id or 'no pose'}" for r in selection.head(8).itertuples()],
        width=960,
        height=540,
    )
