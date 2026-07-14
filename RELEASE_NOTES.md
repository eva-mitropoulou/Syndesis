# Syndesis v1.0.0-paper

This release is the reproducibility package for the Journal of Cheminformatics manuscript, *Pose-coupled native-interaction weighting for kinase ensemble docking: retrospective evaluation on EGFR and CDK2*.

## Scientific scope

- Five EGFR receptor conformations are used for docking.
- The primary EGFR interaction prior uses only four ATP-site holo complexes: 1M17, 1XKK, 4HJO and 5CAV.
- The allosteric ligand JBJ from 6DUK is excluded from ATP-site prior construction; 6DUK remains a receptor conformation.
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
- Replicate-level MD metrics, gates and final decisions.
- Publication figures, Quarto source and rendered PDF.

The release contains computational ranking and pose-persistence evidence. It does not claim experimental activity, affinity, selectivity or current vendor availability.
