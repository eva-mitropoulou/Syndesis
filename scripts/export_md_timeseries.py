"""Export compact per-frame MD evidence from completed production trajectories.

The publication package excludes multi-gigabyte solvated trajectories. This script
exports the exact geometric and key-interaction measurements used by the MD gate
as per-replicate Parquet tables, using the PBC-corrected trajectories generated
during Stage 11 analysis.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import MDAnalysis as mda
import numpy as np
import pandas as pd
from MDAnalysis.analysis import align, rms
from MDAnalysis.lib.distances import distance_array

from syndesis.stage11.interaction_persistence import (
    _interaction_cutoff,
    _select_residue,
)


ROOT = Path(__file__).resolve().parents[1]


def _ligand(universe: mda.Universe):
    selection = universe.select_atoms("resname UNL and not name H*")
    if len(selection) == 0:
        selection = universe.select_atoms("not protein and not resname SOL NA CL and not name H*")
    if len(selection) == 0:
        raise ValueError("ligand atom selection is empty")
    return selection


def export_replicate(metric: dict, interactions: pd.DataFrame, destination: Path) -> tuple[int, int]:
    trajectory = Path(metric["trajectory_file"])
    fitted = trajectory.parent / f"{trajectory.stem}_pbcfit.xtc"
    topology = trajectory.parent / f"{trajectory.stem}_noPBCwater.gro"
    if not fitted.exists() or not topology.exists():
        raise FileNotFoundError(f"missing PBC-corrected inputs for {trajectory}")

    universe = mda.Universe(str(topology), str(fitted))
    protein = universe.select_atoms("protein")
    backbone = universe.select_atoms("protein and backbone")
    ligand = _ligand(universe)
    if len(protein) == 0 or len(backbone) == 0:
        raise ValueError("protein or backbone selection is empty")

    universe.trajectory[0]
    ligand_ref = ligand.positions.copy()
    backbone_ref = backbone.positions.copy()
    reference_pair = distance_array(protein.positions, ligand_ref, box=universe.dimensions)
    pocket_residues = protein[reference_pair.min(axis=1) <= 6.0].residues
    pocket_backbone = pocket_residues.atoms.select_atoms("backbone")
    if len(pocket_backbone) < 4:
        pocket_backbone = backbone
    pocket_ref = pocket_backbone.positions.copy()
    pocket_ref_centered = pocket_ref - pocket_ref.mean(axis=0)
    pocket_atoms = pocket_residues.atoms.select_atoms("not name H*")
    ligand_com_ref = ligand.center_of_mass()

    prepared = []
    for row in interactions.to_dict("records"):
        residue, warnings = _select_residue(universe, row)
        atoms = None if residue is None else residue.atoms.select_atoms("not name H*")
        prepared.append((row, atoms, _interaction_cutoff(str(row["interaction_type"])), warnings))

    geometric_rows = []
    interaction_rows = []
    for frame_index, ts in enumerate(universe.trajectory):
        mobile = pocket_backbone.positions
        mobile_center = mobile.mean(axis=0)
        rotation, _ = align.rotation_matrix(mobile - mobile_center, pocket_ref_centered)
        aligned_pocket = (mobile - mobile_center) @ rotation.T
        aligned_ligand = (ligand.positions - mobile_center) @ rotation.T
        ligand_ref_aligned = ligand_ref - pocket_ref.mean(axis=0)
        in_pocket = True
        if len(pocket_atoms):
            in_pocket = bool(distance_array(pocket_atoms.positions, ligand.positions, box=ts.dimensions).min() <= 4.5)
        common = {
            "md_candidate_id": metric["md_candidate_id"],
            "md_system_id": metric["md_system_id"],
            "replicate_id": metric["replicate_id"],
            "frame_index": frame_index,
            "time_ps": float(ts.time),
        }
        geometric_rows.append({
            **common,
            "ligand_rmsd_angstrom": float(np.sqrt(np.mean(np.sum((aligned_ligand - ligand_ref_aligned) ** 2, axis=1)))),
            "pocket_backbone_rmsd_angstrom": float(np.sqrt(np.mean(np.sum((aligned_pocket - pocket_ref_centered) ** 2, axis=1)))),
            "protein_backbone_rmsd_angstrom": float(rms.rmsd(backbone.positions, backbone_ref, center=False, superposition=False)),
            "ligand_com_drift_angstrom": float(np.linalg.norm(ligand.center_of_mass() - ligand_com_ref)),
            "ligand_radius_gyration_angstrom": float(ligand.radius_of_gyration()),
            "inside_pocket": in_pocket,
        })
        for key, atoms, cutoff, warnings in prepared:
            if atoms is None or len(atoms) == 0:
                minimum = np.nan
                present = False
            else:
                minimum = float(distance_array(atoms.positions, ligand.positions, box=ts.dimensions).min())
                present = bool(minimum <= cutoff)
            interaction_rows.append({
                **common,
                "key_interaction_id": key["key_interaction_id"],
                "residue_role": key["residue_role"],
                "interaction_type": key["interaction_type"],
                "cutoff_angstrom": cutoff,
                "minimum_distance_angstrom": minimum,
                "present": present,
                "selection_warnings": ";".join(warnings),
            })

    geometric_dir = destination / "geometric"
    interaction_dir = destination / "interactions"
    geometric_dir.mkdir(parents=True, exist_ok=True)
    interaction_dir.mkdir(parents=True, exist_ok=True)
    stem = f"{metric['md_system_id']}__{metric['replicate_id']}"
    pd.DataFrame(geometric_rows).to_parquet(geometric_dir / f"{stem}.parquet", index=False)
    pd.DataFrame(interaction_rows).to_parquet(interaction_dir / f"{stem}.parquet", index=False)
    return len(geometric_rows), len(interaction_rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--md-root", type=Path, required=True, help="Directory containing completed Stage 11 trajectories")
    parser.add_argument("--output", type=Path, default=ROOT / "results" / "md" / "timeseries")
    args = parser.parse_args()

    metrics = pd.read_csv(ROOT / "results" / "md" / "md_metrics.csv")
    keys = pd.read_csv(ROOT / "results" / "md" / "md_interaction_persistence.csv")
    manifest_rows = []
    for metric in metrics[metrics["trajectory_analysis_status"].eq("complete")].to_dict("records"):
        trajectory_name = Path(metric["trajectory_file"]).name
        system_dir = args.md_root / metric["md_system_id"]
        metric["trajectory_file"] = str(system_dir / trajectory_name)
        replicate_keys = keys[
            (keys["md_system_id"].eq(metric["md_system_id"]))
            & (keys["replicate_id"].eq(metric["replicate_id"]))
            & (keys["persistence_status"].eq("computed"))
        ].drop_duplicates("key_interaction_id")
        geometric_frames, interaction_rows = export_replicate(metric, replicate_keys, args.output)
        manifest_rows.append({
            "md_system_id": metric["md_system_id"],
            "replicate_id": metric["replicate_id"],
            "geometric_frames": geometric_frames,
            "interaction_rows": interaction_rows,
        })
    pd.DataFrame(manifest_rows).to_csv(args.output / "manifest.csv", index=False)
    print(args.output / "manifest.csv")


if __name__ == "__main__":
    main()
