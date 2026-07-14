"""Audit heavy-atom consistency between docking and ProLIF receptors.

The interaction fingerprint must not use protein heavy atoms that were absent from
the receptor supplied to docking. This script compares the clean docking PDB,
its Open Babel PDBQT derivative, and the PDBFixer hydrogenated ProLIF input for
each selected receptor. It writes a machine-readable, fail-closed audit table.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]


def _atom_key(line: str) -> tuple[str, int, str, str, str]:
    return (
        line[21:22].strip(),
        int(line[22:26].strip()),
        line[26:27].strip(),
        line[17:20].strip(),
        line[12:16].strip(),
    )


def _coordinates(line: str) -> tuple[float, float, float]:
    return (float(line[30:38]), float(line[38:46]), float(line[46:54]))


def _heavy_atoms(path: Path) -> dict[tuple[str, int, str, str, str], tuple[float, float, float]]:
    atoms: dict[tuple[str, int, str, str, str], tuple[float, float, float]] = {}
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if not line.startswith(("ATOM", "HETATM")):
            continue
        element = line[76:78].strip().upper()
        name = line[12:16].strip().upper()
        if element == "H" or (not element and name.startswith("H")):
            continue
        atoms[_atom_key(line)] = _coordinates(line)
    if not atoms:
        raise RuntimeError(f"No protein heavy atoms parsed from {path}")
    return atoms


def _max_coordinate_delta(
    left: dict[tuple[str, int, str, str, str], tuple[float, float, float]],
    right: dict[tuple[str, int, str, str, str], tuple[float, float, float]],
) -> float | None:
    common = left.keys() & right.keys()
    if not common:
        return None
    return max(
        sum((a - b) ** 2 for a, b in zip(left[key], right[key], strict=True)) ** 0.5
        for key in common
    )


def audit(receptor_table: Path, docking_dir: Path, prolif_dir: Path) -> pd.DataFrame:
    receptors = pd.read_parquet(receptor_table)
    rows: list[dict[str, object]] = []
    for record in receptors.to_dict("records"):
        receptor_id = str(record["receptor_id"])
        docking_pdb = docking_dir / f"{receptor_id}.pdb"
        docking_pdbqt = docking_dir / f"{receptor_id}.pdbqt"
        prolif_pdb = prolif_dir / f"docking_receptors__{receptor_id}.h.pdb"
        required = [docking_pdb, docking_pdbqt, prolif_pdb]
        if not all(path.exists() for path in required):
            missing = [str(path) for path in required if not path.exists()]
            raise FileNotFoundError(f"Missing receptor representation for {receptor_id}: {missing}")
        clean = _heavy_atoms(docking_pdb)
        docked = _heavy_atoms(docking_pdbqt)
        prolif = _heavy_atoms(prolif_pdb)
        missing_from_docking = sorted(prolif.keys() - docked.keys())
        missing_from_prolif = sorted(docked.keys() - prolif.keys())
        status = "pass" if not missing_from_docking and not missing_from_prolif else "fail"
        rows.append(
            {
                "receptor_id": receptor_id,
                "docking_pdb": str(docking_pdb),
                "docking_pdbqt": str(docking_pdbqt),
                "prolif_pdb": str(prolif_pdb),
                "clean_pdb_heavy_atoms": len(clean),
                "docking_pdbqt_heavy_atoms": len(docked),
                "prolif_heavy_atoms": len(prolif),
                "prolif_atoms_missing_from_docking_count": len(missing_from_docking),
                "docking_atoms_missing_from_prolif_count": len(missing_from_prolif),
                "clean_to_docking_max_coordinate_delta_angstrom": _max_coordinate_delta(clean, docked),
                "clean_to_prolif_max_coordinate_delta_angstrom": _max_coordinate_delta(clean, prolif),
                "docking_to_prolif_max_coordinate_delta_angstrom": _max_coordinate_delta(docked, prolif),
                "prolif_atoms_missing_from_docking_json": str(missing_from_docking),
                "docking_atoms_missing_from_prolif_json": str(missing_from_prolif),
                "consistency_status": status,
            }
        )
    frame = pd.DataFrame(rows)
    if not frame["consistency_status"].eq("pass").all():
        failed = frame.loc[frame["consistency_status"].ne("pass"), "receptor_id"].tolist()
        raise RuntimeError(f"Receptor heavy-atom consistency audit failed: {failed}")
    return frame


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--receptor-table", type=Path, default=ROOT / "data/processed/stage3/receptor_docking_prep.parquet")
    parser.add_argument("--docking-dir", type=Path, default=ROOT / "data/processed/stage3/docking_receptors")
    parser.add_argument("--prolif-dir", type=Path, default=ROOT / "data/processed/stage5/prolif_receptors")
    parser.add_argument("--output", type=Path, default=ROOT / "results_showcase/submission_robustness/receptor_consistency_audit.csv")
    args = parser.parse_args()
    result = audit(args.receptor_table, args.docking_dir, args.prolif_dir)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(args.output, index=False)
    result.to_parquet(args.output.with_suffix(".parquet"), index=False)
    print(f"wrote {args.output} ({len(result)} receptors; all passed)")


if __name__ == "__main__":
    main()
