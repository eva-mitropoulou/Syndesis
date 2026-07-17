# Frozen results

This directory contains the compact result tables used in the Syndesis paper.
They support inspection and value verification without redistributing raw
docking campaigns or molecular-dynamics trajectories.

- `statistics/paper_metrics.csv` and `statistics/paper_effects.csv` contain
  the primary EGFR and CDK2 ranking values and paired bootstrap effects.
- `statistics/cdk2/` contains the corrected four-receptor CDK2 result tables,
  including the 1PXN preparation summary and late-fusion comparison.
- `cdk2_1pxn_corrected/` contains the corrected 1PXN receptor, selected poses,
  GNINA scores, and interaction fingerprints behind the CDK2 correction.
- `md/persistence_summary.csv` contains the compact system-level summary used
  for the short-timescale modeled-pose persistence analysis.
- `../analysis/pose_coupling_traceability/` contains the frozen tie-aware
  late-fusion traceability audit, its per-ligand and subset tables, and the
  deterministic EGFR representative-case poses and interaction data.

The package is a frozen paper-analysis record. It is not a substitute for the
third-party structural databases, DUD-E source files, docking executables, or
raw MD trajectories required to repeat all upstream calculations from scratch.
