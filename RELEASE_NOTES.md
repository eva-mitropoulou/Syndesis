# Syndesis v1.1.4-paper

This release is the reproducibility package for the Journal of Cheminformatics manuscript, *Pose-coupled native-interaction weighting for kinase ensemble docking: retrospective evaluation on EGFR and CDK2*.

This release standardizes the project, Python package, command-line interface, and
environment name as `Syndesis`. The EGFR primary case study uses the four ATP-site
receptor states; the original 6DUK-inclusive result is retained only as an explicitly
labelled ensemble sensitivity analysis.

This update finalizes the focused manuscript and its reproducibility package. It adds
the DUD-E benchmark limitation, formula-weight sensitivity, CDK2 permutation results,
explicit docking-derived ProLIF chemical perception, table captions, and the corrected
leave-one-receptor-out range. It retains the fail-closed docking/ProLIF receptor-consistency
audit and documents the unphosphorylated 1QMZ representation.

## Scientific scope

- Four EGFR receptor conformations are used for the primary docking and ranking analysis: 1M17, 1XKK, 4HJO, and 5CAV.
- The primary EGFR interaction prior uses the same four ATP-site holo complexes: 1M17, 1XKK, 4HJO and 5CAV.
- Ligand-stripped 6DUK is excluded from the primary protocol and retained only for the five-receptor ensemble sensitivity analysis.
- The primary score couples GNINA CNNscore and interaction recall from the same Uni-Dock-selected pose.
- EGFR is method development plus retrospective evaluation; CDK2 is a target-transfer boundary analysis.

## Archived evidence

- Pose-level EGFR and CDK2 benchmark fingerprints and scores.
- Full 1,000-draw permutation distributions and 2,000-bootstrap summaries.
- Receptor, native-complex, distinct-ligand, exact-overlap and formula sensitivities.
- Native-ligand similarity strata and graph-mapping validation.
- The deterministic 2,000-compound prospective input and corrected ranking.
- Deterministic analog lineage.
- Path-independent GAFF2 ligand parameters, GROMACS topologies and MDP files for all seven reported MD systems.
- Replicate-level MD metrics, gates and final decisions, plus 42 per-frame geometric and interaction time-series tables (42,021 trajectory frames and 588,294 interaction rows).
- Publication figures, Quarto source and rendered PDF.

The release contains computational ranking and pose-persistence evidence. It does not claim experimental activity, affinity, selectivity or current vendor availability.
