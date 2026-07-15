# Frozen benchmark tables

`benchmark/` contains the pose-level tables used to recompute ranking metrics
from the released analysis package. The EGFR and CDK2 tables record benchmark
labels, receptor-specific GNINA scores, native-interaction fingerprints, and
the molecular descriptors used by the matched permutation control.

The input structures and benchmark compounds originate from the Protein Data
Bank and DUD-E. This repository distributes derived research tables for the
paper analysis, not replacement copies of those source resources.

`md_system_code_map.csv` maps the concise MD labels used in the paper (C001–
C003, A004–A006, and N002) to the machine-readable candidate identifiers.
