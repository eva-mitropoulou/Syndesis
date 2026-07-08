# Stage 2 Receptor Ensemble Protocol

Stage 2 selects a compact, state-aware EGFR receptor ensemble from the curated Stage 1
co-crystal benchmark. It does not run docking, GNINA, ProLIF, MD, or analog generation.

The workflow is:

```text
Stage 1 included complexes
-> hard receptor filters
-> pocket residue mapping
-> state-label normalization
-> pocket geometry features
-> state-stratified clustering
-> medoid selection
-> receptor ensemble export
```

Selection is not based on resolution alone. Receptors must pass scope filters, preserve a
usable ATP-site ligand/reference pose, and add pocket or receptor-state diversity.

Stage 3 will validate or prune this ensemble through redocking and cross-docking.

