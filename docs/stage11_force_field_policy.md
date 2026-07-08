# Stage 11 Force-Field Policy

The primary supported stack is CHARMM36m protein plus CGenFF ligand parameters with CHARMM TIP3P water and compatible ions.

The workflow must not mix random force-field families. If CGenFF ligand parameterization is unavailable, setup fails explicitly instead of substituting another ligand force field.

## Amendment (open AMBER/GAFF2 route)

The pipeline now supports and, by default, uses an open force-field stack: **AMBER19SB** (protein) + **GAFF2** (ligand) + **AM1-BCC** (ligand charges) + **OPC3** (water) with AMBER-compatible ions, parameterized via **ACPYPE** (AmberTools: antechamber/parmchk2). This is a deliberate, documented choice made because of CGenFF/ParamChem licensing friction, which makes the CHARMM36m + CGenFF route impractical to run automatically.

This amendment supersedes the "primary supported stack" statement above for default runs: the open AMBER/GAFF2 route is a fully supported alternative rather than a prohibited substitution. The workflow still must not mix force-field families within a single system; the chosen stack (whether CHARMM36m/CGenFF or AMBER19SB/GAFF2) is applied consistently. The force-field choice used for a given system is recorded per-run in the force-field policy output.
