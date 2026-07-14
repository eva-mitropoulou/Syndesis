# Analysis-ready results

This directory contains the compact, machine-readable results used by the
manuscript. It intentionally excludes raw docking poses and MD trajectories.

## `analysis_inputs/`

- `egfr_master.parquet` and `cdk2_master.parquet`: one strict pose-level row per benchmark molecule and receptor state.
- `*_native_interaction_fingerprints.parquet`: native-complex fingerprints used to define each prior.
- `*_ligand_descriptors.csv`: benchmark identifiers, labels, heavy-atom counts, and molecular weights used by matched nulls.
- `prospective_ranked_corrected.parquet`: selected-pose prospective ranking used for the gate audit and figure.
- `permutation_null_draws.parquet`: archived permutation draws used in the density figure.

These derived tables permit the statistical analyses and figures to be rerun
without redistributing third-party benchmark structures or raw trajectories.
Machine-specific pose paths and execution-log payloads are excluded; identifiers,
scores, labels, fingerprints, interaction counts, and status fields are retained.

## `robustness/`

- `bootstrap_metric_intervals.csv`: class-stratified 95% intervals for each target and ranking arm.
- `paired_metric_effects.csv`: paired coupled-minus-GNINA effects.
- `permutation_null_summary.csv`: 1,000-draw unrestricted, property-matched, and class-conditional nulls.
- `leave_one_receptor_out.csv`: receptor-ensemble sensitivity.
- `*_native_prior_sensitivity.csv`: prior-definition and native-ligand exclusion analyses.
- `*_native_interaction_bits.csv`: exact residue-by-interaction-type prior definitions.
- `interaction_formula_sensitivity.csv`: interaction term, combination rule, and weight ablations.
- `interaction_size_correlations.csv`: recall correlations with molecular size and contact count.
- `top1_active_counts.csv`: raw active counts at the 1% cutoff.
- `prospective_gate_audit.csv`: prospective thresholds and selected-pose gate audit.
- `prospective_gate_pass_candidates.csv`: gate-passing, unlabeled screening hypotheses.
- `deterministic_analog_lineage.csv`: parent, transformation, and pose provenance.

## `md/`

- `ligand_parameterization_report.csv`: GAFF2/AM1-BCC parameterization status and warnings.
- `system_builds.csv` and `md_runs.csv`: system composition and replicate execution records.
- `md_metrics.csv`, `md_interaction_persistence.csv`, and
  `md_binding_mode_persistence.csv`: replicate-level trajectory measurements.
- `md_candidate_summary.csv`: candidate-level majority-replicate verdicts.

All rankings are computational prioritizations. Prospective rows are not labeled
active, and MD stability is not interpreted as binding affinity.
