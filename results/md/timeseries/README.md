# Frame-Level MD Evidence

This directory contains the compact per-frame measurements used to derive the
published MD summaries and stability decisions. Each production replicate has:

- `geometric/<system>__<replicate>.parquet`: ligand RMSD in the pocket-aligned
  frame, pocket and protein RMSD, ligand center-of-mass drift, radius of
  gyration, and pocket-retention status.
- `interactions/<system>__<replicate>.parquet`: minimum ligand--residue distance
  and present/absent status for every key interaction at every saved frame.

`manifest.csv` records the frame and interaction-row counts. The tables were
exported from the completed PBC-corrected, water-stripped Stage 11 trajectories.
They reproduce the median RMSD and interaction occupancy values in the MD
summary tables exactly. The original solvated trajectories are excluded because
they are large generated artifacts; the topology, parameters, compact input
coordinates, and deterministic export script are retained in this release.
