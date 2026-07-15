# Syndesis v1.2.0-paper

This is the frozen reproducibility package for *Pose-coupled native-interaction
weighting in kinase ensemble docking: EGFR enrichment and CDK2 transfer*.

## Primary analyses

- EGFR uses four ATP-site receptors and four matching native complexes: 1M17,
  1XKK, 4HJO, and 5CAV. Ligand-stripped 6DUK is retained only for an
  ensemble-sensitivity analysis.
- CDK2 uses four docking receptors and four matching native complexes: 1FIN,
  2A4L, 1AQ1, and 1PXN. 1QMZ is excluded from both primary components.
- The corrected 1PXN preparation retains Lys33-A (occupancy 0.98), removes
  Lys33-B (occupancy 0.02), and has identical docking and ProLIF heavy-atom
  representations.
- CDK2 EF1% increases from 11.60 for GNINA to 14.98 for pose coupling; the
  paired difference is 3.375 (95% CI 0.422–5.695).

## Included evidence

- Corrected 1PXN receptor, top poses, GNINA scores, strict ProLIF
  fingerprints, and receptor-consistency audit.
- EGFR and CDK2 score tables, bootstrap and permutation outputs, and
  receptor/native-prior sensitivity results.
- Publication figures, Quarto manuscript source, rendered PDF, and
  machine-readable MD provenance and gate outputs.

Large campaign intermediates, trajectories, and unrelated prospective-screening
materials are excluded. This release reports computational ranking and
short-timescale pose-persistence evidence; it does not claim experimental
activity, affinity, selectivity, or clinical relevance.
