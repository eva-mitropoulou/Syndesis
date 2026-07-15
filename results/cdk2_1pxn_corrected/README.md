# Corrected 1PXN rerun outputs

This directory contains the targeted 1PXN rerun used by the primary four-receptor,
four-native-complex CDK2 analysis. The receptor retains the 0.98-occupancy Lys33-A
side chain and excludes the 0.02-occupancy Lys33-B conformer.

- `1pxn_a_ck6_corrected_poses.tar.gz`: corrected receptor PDBQT, GNINA pose SDF,
  and one Uni-Dock-selected top pose per ligand.
- `scores_1pxn_a_ck6.parquet`: GNINA score-only outputs for the frozen 28,296
  CDK2 molecules.
- `pose_fingerprints_1pxn_a_ck6.parquet`: strict graph-preserving ProLIF
  fingerprints for the same poses.
- `1pxn_a_lys33a_prolif_receptor.pdb`: the ProLIF receptor representation.

The corresponding four-receptor score table, paired bootstrap, permutation draws,
and sensitivity outputs are in `../robustness/cdk2_four_receptor_four_prior/`.
