# Stage 11 Ligand Parameterization

Ligands use the exact protonation and tautomer states selected by previous stages. CGenFF `.str` files may be supplied under `data/inputs/stage11/cgenff_str/`.

Missing `.str` files, missing CGenFF, or missing `cgenff_charmm2gmx` conversion are recorded as parameterization failures. High CGenFF penalty scores must be flagged before flagship MD.
