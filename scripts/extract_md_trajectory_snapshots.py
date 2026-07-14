#!/usr/bin/env python3
"""Export the first and last saved frames of a selected MD production run."""

from pathlib import Path

import MDAnalysis as mda


OUTPUT_DIR = Path("figures/manuscript/md_trajectory_frames")

SYSTEMS = (
    ("control_002_rep01", Path("/mnt/e/egfr_md_work/mdsys_mdcand_002")),
    ("misdocked_control_rep01", Path("/mnt/e/egfr_md_work/mdsys_mdcand_neg01")),
)


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    rows = ["label,frame,time_ps,coordinate_file"]
    for system_label, system_root in SYSTEMS:
        universe = mda.Universe(
            str(system_root / "production_quick_noPBCwater.gro"),
            str(system_root / "production_quick_pbcfit.xtc"),
        )
        selection = universe.select_atoms("protein or resname UNL")
        snapshots = (("start", 0), ("end", len(universe.trajectory) - 1))
        for label, frame in snapshots:
            timestep = universe.trajectory[frame]
            output = OUTPUT_DIR / f"{system_label}_{label}.pdb"
            selection.write(str(output))
            rows.append(
                f"{system_label}_{label},{frame},{timestep.time:.3f},{output.as_posix()}"
            )

    (OUTPUT_DIR / "control_002_rep01_snapshot_manifest.csv").write_text(
        "\n".join(rows) + "\n", encoding="utf-8"
    )


if __name__ == "__main__":
    main()
