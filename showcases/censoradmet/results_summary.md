# CensorADMET Results Summary

## Research question

ADMET assays often report detection limits and inequality labels. Treating a
bound as an exact value can create systematic errors, while discarding censored
rows throws away information. CensorADMET evaluates these choices under both
ordinary and chemically shifted test distributions.

## Benchmark construction

The released benchmark contains 31 ChEMBL-derived endpoint datasets with
39,157 curated parent-compound records. It preserves 21,922 exact labels and
17,235 one-sided or interval-censored labels. Eight endpoints satisfy the
prespecified sample, censoring, structure, and metadata requirements for the
main modelling analysis. The central modelling subset is dominated by CYP and
hERG endpoints; conclusions for permeability, solubility, clearance, and
transport remain open questions.

## Evaluation design

The analysis compares two common naive treatments of censoring with interval-
aware approaches:

- threshold-as-exact XGBoost;
- dropping censored observations;
- XGBoost accelerated-failure-time regression;
- a censored Gaussian Tobit ensemble.

Models use Morgan fingerprints for the main paired comparisons. Evaluation uses
random, scaffold, Butina-cluster, and publication-year pseudo-temporal splits.
Metrics include exact-label MAE, Spearman correlation, one-sided bound violation,
calibration coverage, interval width, and natural-label bound compatibility.

## Main result

The endpoint-level analysis uses eight analysis units rather than treating
correlated seeds and splits as independent observations. Relative to naive
baselines, interval-aware methods reduce one-sided violations by 0.61--0.78,
but increase exact-label MAE by 0.68--0.85 p-units and reduce Spearman
correlation by 0.12--0.17. The direction is consistent across the eight
endpoints, with false-discovery-rate-adjusted paired-test values of 0.0078.

This is a practical decision result: a workflow that must obey assay bounds may
prefer interval-aware objectives, while a workflow prioritising exact-label
ranking may prefer a simpler baseline.

## Calibration result

At 90% nominal coverage, absolute conformal calibration reaches 0.908 empirical
exact-label coverage with a mean interval width of 5.111 p-units. The width is
large: five p-units span five orders of magnitude on the concentration scale.
Calibration therefore recovers empirical coverage while also demonstrating
that uncertainty estimates may be too broad for some medicinal-chemistry
decisions. Natural censored labels are evaluated by bound compatibility, not by
claiming hidden true-value coverage.

## Why this is technically interesting

The project demonstrates an end-to-end cheminformatics workflow:

- structured extraction and curation of heterogeneous bioactivity records;
- explicit interval-label semantics instead of silent threshold conversion;
- chemical-distribution shift evaluation beyond a single random split;
- numerically stable censored Gaussian likelihoods evaluated in log space;
- endpoint-level statistical inference with bootstrap confidence intervals;
- calibration analysis that separates exact, synthetic hidden-truth, and natural
  censored-label evidence;
- release-level numerical checks connecting source metrics, tables, and figures.

The central contribution is diagnostic: it makes the accuracy-versus-constraint
trade-off visible and reproducible.
