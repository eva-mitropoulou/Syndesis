# Syndesis

**Pose-coupled native-interaction weighting for auditable kinase docking**

[![Python 3.11+](https://img.shields.io/badge/python-3.11%2B-3776AB.svg)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/license-MIT-2F855A.svg)](LICENSE)

Syndesis is a structure-based cheminformatics workflow for testing whether interactions observed in native protein-ligand complexes add useful ranking information to neural docking scores. Its central rule is deliberately simple: combine GNINA CNNscore with recovery of a target-native interaction prior **on the same receptor-specific pose**, then maximize that coupled score across an ensemble.

I built Syndesis to make virtual-screening decisions inspectable. A high score is not accepted as a biological claim; the workflow records the pose, recovered contacts, uncertainty, native-ligand overlap, parameterization warnings, and replicate-level molecular-dynamics evidence behind each decision. The name comes from the Greek *syndesis*, or “connection,” reflecting both molecular binding and the connection of independent evidence.

> **Claim boundary:** this repository reports computational ranking and pose-persistence evidence. It does not claim experimental activity, affinity, selectivity, or clinical relevance.

## Main Result

On the EGFR DUD-E benchmark (542 actives and 35,010 decoys), strict graph-preserving ProLIF recomputation produced:

| Ranking | ROC-AUC | EF1% | EF5% | BEDROC |
|---|---:|---:|---:|---:|
| GNINA | 0.766 | 11.79 | 6.90 | 0.208 |
| GNINA + ATP-site native-union recall | **0.772** | **16.40** | **7.75** | **0.282** |

At the 1% cutoff, the coupled ranking retrieved **89 actives among 356 molecules**, compared with **64** for GNINA. The paired EF1% gain was 4.61 (95% bootstrap CI 2.58–6.82).

The EGFR docking ensemble contains five receptor conformations, including 6DUK. The primary interaction prior contains only the four ATP-site holo complexes 1M17, 1XKK, 4HJO, and 5CAV. JBJ from 6DUK is allosteric and is excluded from every ATP-site prior calculation.

The result was tested against three different 1,000-permutation nulls:

| Null | Mean EF1% | Observed minus null | Empirical p |
|---|---:|---:|---:|
| All-ligand shuffle | 11.35 | 5.05 | 0.0010 |
| Heavy-atom-count-matched shuffle | 12.39 | 4.01 | 0.0010 |
| Class-conditional assignment | 14.26 | 2.13 | 0.0040 |

Every leave-one-receptor-out EGFR analysis preserved a positive paired effect. Excluding the exact-overlap FMM complex retained EF1% 15.29, and removing both duplicate AQ4 complexes together retained EF1% 15.66; both paired intervals excluded zero. Among the 369 actives with ECFP4 similarity below 0.30 to every distinct ATP-site native ligand, the coupled score recovered 54 in the global top 1%, compared with 39 for GNINA.

CDK2 defined the boundary of the claim. The same fixed rule increased EF1% from
10.97 to 13.08, retrieving 62 rather than 52 actives among the first 283
molecules. Its paired EF1% difference was 2.11 (95% CI -0.42 to 4.64), and
leave-one-receptor-out effects were heterogeneous; the result is therefore
reported as favorable but unresolved rather than as an independent replication.

## What Is Different

The primary score for ligand `i` and receptor state `r` is:

```text
S(i) = max_r CNNscore(i,r) * [1 + recall(i,r)]
```

where `recall(i,r)` is the fraction of target-native residue-by-interaction bits recovered by that exact pose. CNNscore from one receptor is never combined with interaction recall from another.

The implementation is fail-closed:

- Docked coordinates are mapped back onto the prepared SDF molecular graph before ProLIF analysis.
- Atom-count, element-order, or fingerprint failures block ranking instead of becoming zero-valued interactions.
- Missing ligand net charge blocks MD parameterization.
- GAFF2 ligands cannot be mixed with a CHARMM protein force field unless explicitly overridden.
- Parameterization warnings are retained in machine-readable reports.
- No dependency or scientific-method fallback is selected silently.

## Evidence Layers

1. **Structural curation:** ATP-site kinase complexes, normalized residue maps, and receptor-state metadata.
2. **Ensemble docking:** five EGFR or CDK2 receptor conformations with reproducible seeds and boxes.
3. **Neural rescoring:** GNINA score-only evaluation of one pose per ligand-receptor pair.
4. **Interaction analysis:** ProLIF fingerprints using prepared ligand graphs and docked coordinates.
5. **Statistical controls:** paired bootstraps, three permutation nulls, formula sensitivity, and complete leave-one-out analyses.
6. **Leakage audit:** exact identity, scaffold, and ECFP4 similarity between native ligands and benchmark actives.
7. **Deterministic analog audit:** versioned RDKit transformations with complete parent-to-analog lineage.
8. **Replicated MD:** 7 systems x 3 independent 20 ns trajectories using ff19SB, GAFF2/AM1-BCC, OPC3, and GROMACS, with frame-level geometric and interaction time series.

The MD gate accepted four systems by majority replicate and rejected three, including a deliberately mis-docked negative control with median ligand RMSD 5.72 Å and key-interaction occupancy 0.008.

## Repository Map

```text
src/syndesis/   installable pipeline and CLI
configs/                 versioned scientific policies and run parameters
scripts/                 enrichment, robustness, provenance, and figure workflows
tests/                   unit and regression tests
docs/                    methods, source registries, and non-claim policies
results/                 compact machine-readable publication results
figures/                 publication figures
manuscript/              Quarto source and rendered paper
workflows/               Snakemake orchestration
```

Large docking intermediates and trajectories are intentionally excluded. The release contains compact results, source manifests, checksums, and deterministic analysis scripts rather than multi-gigabyte run directories.

## Installation

Create the scientific environment and install the package:

```bash
micromamba create -n syndesis -f environment.yml
micromamba activate syndesis
python -m pip install -e .
pytest -q
```

The default suite excludes tests that require generated campaign artifacts or
external chemistry executables. Run `pytest -m integration` after staging those
inputs and tools.

The full workflow also requires working installations of Uni-Dock, GNINA 1.3.3, Open Babel, AmberTools (`antechamber` and `parmchk2`), ACPYPE, and GROMACS 2026.0. Tool paths are deployment settings; see [`configs/tools.example.yaml`](configs/tools.example.yaml) and [`docs/reproducibility_notes.md`](docs/reproducibility_notes.md).

## Reproducing the Paper

The paper analysis is frozen in [`configs/paper_analysis.yaml`](configs/paper_analysis.yaml). Four reproduction levels are supported:

```bash
# 1. Validate package behavior
pytest -q

# 2. Rebuild the ATP-site prior sensitivities and publication statistics
PYTHONPATH=src python scripts/native_prior_sensitivity_analysis.py
PYTHONPATH=src python scripts/submission_robustness_analysis.py

# 3. Rebuild the prospective ranking and native-chemotype audit
PYTHONPATH=src python scripts/build_corrected_prospective_ranking.py
PYTHONPATH=src python scripts/active_native_similarity_analysis.py

# 3b. Re-export frame-level MD evidence from completed trajectories
PYTHONPATH=src python scripts/export_md_timeseries.py --md-root /path/to/egfr_md_work

# 4. Render the manuscript after regenerating figures
python scripts/create_manuscript_figures.py
quarto render manuscript/syndesis_jcheminformatics_v2.qmd --to typst
```

The complete docking campaigns require the source DUD-E and ZINC files described in the manifests, prepared receptor files, and GPU-capable external tools. Commands do not substitute a different engine when a required tool is missing.

## Citation

Please cite the frozen software and data release until the ChemRxiv preprint is posted:

> Mitropoulou E, Giannopoulos D. **Syndesis v1.1.0-paper: pose-coupled native-interaction weighting for kinase ensemble docking.** GitHub release, 2026.

Machine-readable citation metadata are provided in [`CITATION.cff`](CITATION.cff). The reproducibility package is frozen in the `v1.1.0-paper` GitHub release. ChemRxiv preprint metadata will be added after posting.

## Authors

**Evangelia Mitropoulou** and **Dimitris Giannopoulos**<br>
Department of Chemistry, University of Patras, Greece

Syndesis was developed as a research-grade computational chemistry project with emphasis on reproducibility, failure-aware engineering, and scientifically bounded interpretation.
