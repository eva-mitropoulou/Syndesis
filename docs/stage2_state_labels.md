# Stage 2 State Labels

Stage 2 imports KLIFS and KinCore labels from Stage 1 when available. Missing labels are
kept as null and warnings are recorded. The code may compute simple geometric fallback
features, but it does not invent authoritative kinase-state labels.

Default strata:

- `active_like`: active-like or DFGin active controls.
- `inactive_like`: inactive-like controls such as lapatinib-bound EGFR.
- `dfgout_or_typeII_like`: reserved for DFGout or type-II-like structures if present.
- `unknown_state`: used when state metadata are insufficient.

Mutation status is metadata only and is not used to make mutant-selectivity claims.

