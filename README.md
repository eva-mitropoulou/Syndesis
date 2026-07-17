# Syndesis

**Pose-coupled native-interaction weighting for kinase ensemble docking**

Syndesis is an open, reproducible analysis package for testing whether
target-native protein-ligand interactions add early-ranking information to a
neural docking score. Its central rule is simple: combine GNINA CNNscore and
native-interaction recall for the **same receptor-specific pose**, then retain
the highest coupled score across the receptor ensemble.

The package accompanies the paper *Pose-coupled native-interaction weighting in
kinase ensemble docking: EGFR enrichment and CDK2 transfer.* It reports
computational ranking and short-timescale modeled-pose persistence only; it
does not claim experimental activity, affinity, inhibition, or clinical value.

## Main results

| Target | GNINA EF1% | Pose-coupled EF1% | Paired difference (95% CI) |
|---|---:|---:|---:|
| EGFR | 11.98 | 16.40 | 4.42 (2.58–6.63) |
| CDK2 | 11.60 | 14.98 | 3.375 (0.422–5.695) |

On EGFR, the coupled ranking recovered 89 known actives in the first 356
molecules, compared with 65 for GNINA. EGFR was the method-development target.
The fixed rule was then applied to CDK2, where the favorable result was more
sensitive to receptor and prior composition.

The pose-decoupled late-fusion comparator recovered 86 EGFR actives (EF1%
15.85) and 78 CDK2 actives (EF1% 16.45). It is not structurally equivalent to
same-pose coupling: 56.69% of EGFR and 62.02% of CDK2 late-fusion scores combine
maxima from different receptor-specific poses. Coupling therefore preserves
score provenance rather than universally maximizing enrichment.

The EGFR analysis uses four ATP-site structures (1M17, 1XKK, 4HJO, and 5CAV)
for both docking and the native-interaction reference. The CDK2 analysis uses
1FIN, 2A4L, 1AQ1, and 1PXN. The corrected 1PXN preparation retains the
0.98-occupancy Lys33-A conformer, removes the 0.02-occupancy Lys33-B conformer,
and has identical docking and ProLIF heavy-atom representations.

## Included material

```text
src/syndesis/       score and interaction-prior utilities
scripts/            one verification entry point
data/               frozen pose-level benchmark input tables
results/            frozen statistical and MD result tables
figures/            publication figures
manuscript/         Quarto source, references, figures, and rendered PDF
analysis/            pose-coupling traceability audit and representative data
tests/              unit tests for the interaction-prior calculations
```

The repository intentionally excludes internal project notes, pipeline-stage
artifacts, raw docking campaigns, trajectories, temporary render files, and
intermediate workspace records. The released tables are the compact inputs and outputs
needed to inspect the reported analyses.

## Verify the reported values

Python 3.11+ with NumPy, pandas, and PyArrow is sufficient to verify the
released results.

```bash
python -m pip install -e .
make test
make verify
```

`make verify` reads the frozen result tables and checks the reported EGFR and
CDK2 EF1% values and paired effects. It does not rerun docking, GNINA, or MD.

## Citation

Please cite the ChemRxiv preprint:

> Mitropoulou, E.; Giannopoulos, D. *Pose-coupled native-interaction weighting
> in kinase ensemble docking: EGFR enrichment and CDK2 transfer.* ChemRxiv,
> 2026. https://doi.org/10.26434/chemrxiv.15006204/v1

Citation metadata are provided in [CITATION.cff](CITATION.cff). The immutable
paper snapshot is [v1.3.0-paper](https://github.com/eva-mitropoulou/Syndesis/tree/v1.3.0-paper).

## License

Code and the released derived analysis tables are available under the
[MIT License](LICENSE). The source DUD-E and Protein Data Bank records remain
subject to their respective terms and attribution requirements.
