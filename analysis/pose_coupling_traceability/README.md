# Pose-coupling traceability audit

This directory audits the structural interpretability of the pose-decoupled
late-fusion control. It uses frozen, receptor-level scores and ProLIF
fingerprints only; it does not invoke Uni-Dock, GNINA, ProLIF, molecular
dynamics, bootstrapping, or permutations.

## Inputs

- The compact public package uses `../../data/benchmark/egfr_pose_scores.parquet`,
  `../../data/benchmark/cdk2_pose_scores.parquet`, and
  `../../data/benchmark/egfr_native_fingerprints.parquet`.
- The working manuscript repository uses the corresponding frozen master tables
  and native-fingerprint table under `results_showcase/` and `data/processed/`.

In both layouts, the analysis filters to the four primary EGFR receptors
`1M17`, `1XKK`, `4HJO`, and `5CAV`, and the four primary CDK2 receptors
`1FIN`, `2A4L`, `1AQ1`, and corrected `1PXN`.

Ligands use the stable lexical order created by the frozen manuscript
`pivot()` analyses. Score ties are ranked with stable descending mergesort;
receptor maxima are analysed tie-aware with `numpy.isclose(rtol=1e-10,
atol=1e-12)` and a fixed primary receptor order only for the supplementary,
deterministic tie-break label.

## Run

```bash
python analysis/pose_coupling_traceability/traceability_analysis.py
```

The default run regenerates the numerical audit only. In the full working
repository, pass `--render` to regenerate the selected-pose PyMOL scripts and
images from the source pose files. The tagged public package includes the
frozen PML scripts and selected pose files; no workstation-specific paths are
written to the reusable outputs.

## Outputs

- `traceability_summary.csv`: all-molecule pose-realizability and fusion-gap summaries.
- `traceability_by_ligand.csv`: one row per target–ligand, scores, ranks,
  tie-aware maxima, top-1% membership, coverage, and rescued/lost labels.
- `traceability_subsets.csv`: the requested molecule, active, decoy, and
  ranked-subset counts.
- `fusion_gap_summary.csv`: descriptive late-fusion-minus-coupled gap statistics.
- `maxima_tie_summary.csv`: unique and tied receptor maxima under the tie-aware rule.
- `rescued_lost_shared_actives.csv`, `rescued_actives.csv`, and
  `rescued_active_receptor_transitions.csv`.
- `representative_selection.csv`, `representative_case.csv`, and
  `representative_case_interactions.csv`.
- `representative_case.pml` and, when PyMOL is available,
  `representative_case.png`.

The representative active is selected mechanically: it must be coupled-top-1%,
GNINA-rescued, and pose-non-realizable by late fusion. Cases whose coupled
receptor differs from both independently maximizing receptors are preferred;
the largest GNINA-to-coupled rank improvement then wins, with stable ligand
order resolving any remaining tie.
