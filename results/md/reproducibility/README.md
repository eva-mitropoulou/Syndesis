# MD Reproducibility Inputs

This directory contains the exact compact GROMACS and GAFF2 inputs for the seven systems reported in the paper. Each candidate directory contains:

- the pre-solvation protein-ligand complex coordinates;
- a standalone ff19SB protein topology with a relative ligand-parameter include;
- the GAFF2/AM1-BCC ligand ITP, GRO, MOL2, frcmod and Amber topology;
- protein position restraints;
- minimization, NVT, NPT and three production-replicate MDP files;
- ACPYPE and SQM records used to audit parameter generation.

`checksums.csv` records SHA-256 hashes and file sizes. Replicate-level outcomes are in the parent `results/md/` tables. Solvated coordinates, compiled TPR files and trajectories are excluded because they are large generated artifacts; the system-build table records their composition and the packaged inputs rebuild them deterministically with the documented GROMACS version.
