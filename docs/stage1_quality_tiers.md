# Stage 1 Quality Tiers

## Tier A

Tier A complexes are primary references for redocking and interaction recovery:

- X-ray diffraction preferred.
- Resolution <= 2.5 A.
- Noncovalent ATP-site ligand.
- Ligand heavy atoms >= 12.
- Ligand coordinates complete enough for native-pose reference use.
- Ligand occupancy preferably >= 0.8.
- No severe alternate-conformer ambiguity.
- No missing key active-site residues.
- Validation files available where RCSB provides them.
- Mutation status recorded.
- Receptor-state metadata available or derivable.

## Tier B

Tier B complexes are useful for redocking/cross-docking but have minor caveats:

- Resolution > 2.5 A and <= 3.0 A, or minor local-quality caveats.
- Ligand is usable.
- Active site is mostly complete.
- Ligand remains in scope.

## Tier C

Tier C complexes are retained for state-diversity analysis or fallback use:

- Resolution > 3.0 A and <= 3.5 A, or unique receptor-state value with weaker quality.
- Ligand is usable and in scope.
- Not a primary training reference unless explicitly justified.

## Rejected

Rejected complexes violate scope or hard quality filters. They remain machine-readable in
`data/interim/stage1/rejected_complexes.parquet` with explicit exclusion reasons.

