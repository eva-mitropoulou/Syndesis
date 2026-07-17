# Syndesis v1.3.0-paper

Frozen reproducibility package for *Pose-coupled native-interaction weighting
in kinase ensemble docking: EGFR enrichment and CDK2 transfer*.

- The rendered manuscript and its Quarto source are included under
  `manuscript/`.
- CDK2 uses the four-receptor/four-native-complex primary design: 1FIN, 2A4L,
  1AQ1, and corrected 1PXN. GNINA EF1% is 11.60 and pose-coupled EF1% is
  14.98 (paired difference 3.375; 95% CI 0.422–5.695).
- Pose-decoupled late fusion recovered 78 CDK2 actives at EF1% 16.45. The
  coupled-minus-late-fusion EF1% difference is −1.477 (95% CI −3.797–0.422).
- `analysis/pose_coupling_traceability/` provides the tie-aware audit of the
  late-fusion control, including per-ligand results, fusion-gap summaries,
  representative-case data, selected poses, recovered interaction bits, and
  portable PyMOL scripts.

Large raw docking campaigns and molecular-dynamics trajectories are excluded;
the compact frozen input tables and outputs needed to inspect the reported
analyses are included.
