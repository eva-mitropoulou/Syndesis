# Syndesis

Syndesis is a structure-based computational workflow for EGFR ATP-site analog design that tries to keep only the candidates that bind in a believable way, rather than the ones that simply score well. I built it because a docking score on its own is easy to "hack": a molecule can be scored highly while sitting in a pose that breaks the interactions a real EGFR inhibitor is known to make. Syndesis instead asks several independent lines of structural evidence — the docked pose, a residue-level interaction fingerprint, a leakage-controlled pose-confidence model, and a short molecular-dynamics pose-stability test — to agree before a candidate is accepted. The name (Greek *sýndesis*, "binding together") refers both to ligand–target binding and to binding these independent checks together into one decision.

This is a computational, methods-oriented project. It does **not** claim experimentally confirmed inhibitors, cellular activity, selectivity, or clinical relevance. Scores, poses, and simulations are hypotheses that would require experimental validation.

## Table of Contents

- [What the workflow does](#what-the-workflow-does)
- [Pipeline](#pipeline)
- [Main results](#main-results)
  - [Redocking validation](#redocking-validation)
  - [Pose-confidence model](#pose-confidence-model)
  - [Analog optimization and score-hacking](#analog-optimization-and-score-hacking)
  - [MD pose stability](#md-pose-stability)
- [Scope and limits](#scope-and-limits)
- [Reproduce](#reproduce)
- [Repository layout](#repository-layout)
- [Citation](#citation)

## What the workflow does

The target is the human EGFR (ERBB1/HER1) kinase domain (UniProt P00533), orthosteric ATP site, reversible noncovalent inhibitors. The starting point is a set of curated EGFR co-crystal structures: I extract the native ligands, clean the receptors, and use the crystallographic binding geometry as the reference for what a real ATP-site pose looks like.

From there the idea is simple to state. Most weak docking pipelines stop at "dock, then sort by score". Syndesis goes further: it checks whether a docked pose reproduces the conserved EGFR interactions (hinge Met793, gatekeeper Thr790, catalytic Lys745, DFG Asp855), scores pose confidence with a model that is explicitly audited for label leakage, and — for the finalists — runs short explicit-solvent MD to see whether the pose actually holds. A candidate is accepted only if these agree; a candidate whose score improves while its binding mode degrades is flagged and rejected as score-hacking. The headline metric is therefore not the best docking score but the **accepted-analog rate** under all of these gates.

The project also includes a small, honest experiment on the analog-generation step itself: I compare a deterministic RDKit rule-based generator against a local LLM-agent design loop (Qwen3-32B served through vLLM), screening both through exactly the same acceptance gates.

## Pipeline

The workflow is organised as thirteen configurable stages, each driven from a single CLI (`egfrforge <command> --config configs/<stage>.yaml`) and each writing versioned tables plus an HTML report:

0. Scope definition (target, inclusion/exclusion rules).
1. EGFR co-crystal benchmark: fetch structures, extract native ligands, clean receptors.
2. Receptor ensemble: cluster pocket conformations, select a representative state.
3. Redocking / cross-docking: recover native poses; symmetry-aware RMSD; pose labeling.
4. GNINA rescoring: CNNscore / CNNaffinity as pose-quality and affinity features.
5. Interaction atlas (ProLIF): define the conserved-core EGFR key interactions.
6. Pose-confidence model: predict native-like poses from leakage-audited, deployment-safe features.
7. Candidate library: assemble known EGFR ligands and native ligands.
8. Candidate screening: dock, rescore, interaction-profile, aggregate.
9. Analog optimization: RDKit rule-based and LLM-agent generation, screened head-to-head.
10. Ablation / benchmark: compare strategies by accepted-analog rate and score-hacking rate.
11. MD stress test: short, replicated explicit-solvent MD for finalist pose stability.
12. Candidate dossiers: per-candidate evidence cards, provenance bundle, model and dataset cards.

## Main results

All numbers below come from the pipeline's own outputs. The full replicated finalist MD was still running at the time of writing, so the MD numbers for the three known controls are marked as preliminary and the designed analogs are noted as scheduled.

### Redocking validation

Across the benchmark ligands, docking recovers a near-native pose (≤ 2 Å) among the sampled poses for **5 of 5** ligands, with a median best-pose RMSD of **1.43 Å**. However, the top-scored (rank-1) pose is not always the native-like one — for example the 1M17/AQ4 complex places a 1.6 Å pose at rank 3 while rank 1 is ~8 Å. In other words, docking *samples* the right geometry but the empirical score does not reliably *rank* it first. This is exactly the gap the interaction constraints and the pose-confidence model are meant to close, and it is the reason the workflow does not trust docking score alone.

### Pose-confidence model

The pose-confidence label is the RMSD-to-crystal ground truth (is the pose native-like?), kept deliberately separate from the interaction features so the model cannot read its own answer. Features are passed through a default-deny leakage audit: only deployment-safe inputs (raw interaction bits, docking/GNINA scores, pose geometry, ligand descriptors, receptor-state metadata) are allowed, and any feature that measures agreement with the native reference is dropped.

On a Bemis-Murcko scaffold-holdout test set, the learned ranker (LightGBM LambdaRank) reaches **NDCG@3 = 0.69** versus **0.41** for ranking by docking score alone and 0.28–0.32 for GNINA CNN scores. I treat this as encouraging but not decisive: the benchmark contains only about four distinct scaffolds in the holdout, so I report the difference as a trend and flag the small scaffold diversity as a limitation rather than claiming a powered result.

### Analog optimization and score-hacking

Starting from three known EGFR chemotypes as seeds, the workflow generated and screened 26 analogs and accepted 10 under the full gate set. Screening was identical across generation strategies:

| Strategy | Generator | Valid | Screened | Accepted |
|---|---|---:|---:|---:|
| rule-based | RDKit | 17 | 17 | 6 |
| single agent | Qwen3-32B | 2 | 2 | 1 |
| council (no feedback) | Qwen3-32B | 3 | 3 | 1 |
| council (tool feedback) | Qwen3-32B | 2 | 2 | 1 |
| council (interaction-constrained) | Qwen3-32B | 2 | 2 | 1 |

The LLM arms each produced only a handful of molecules, so I present this as a **feasibility demonstration** — a local open-weights agent can be bridged to the same cheminformatics tooling and produce valid, screenable analogs evaluated on equal footing with the rule-based baseline — not as a powered comparison between generators.

Score-hacking is reported with two complementary definitions to be transparent: a strict acceptance gate flagged and rejected **3** analogs whose CNN score improved while the binding mode broke, and a broader diagnostic auditor flagged **8** cases in total (the 3 plus 5 milder cases where a composite score improved while interaction evidence weakened — 5 of which still passed the acceptance gate). Reporting both is the honest framing: the auditor is stricter than the gate and surfaces borderline candidates the gate lets through.

### MD pose stability

Finalists (three known controls plus the top three accepted analogs) are stress-tested with short explicit-solvent MD (AMBER19SB / GAFF2 / AM1-BCC / OPC3, GROMACS), with restrained equilibration and three independent replicates per finalist. Pose stability is measured in the protein reference frame after periodic-boundary correction and least-squares fitting on the protein backbone; ligand RMSD is taken relative to the docked start pose, with interaction persistence tracked alongside.

Preliminary single-replicate results for the three known EGFR chemotypes show stable poses over 20 ns (median ligand RMSD ≈ 1.8–2.0 Å, pocket contact retained throughout). A methodological note worth flagging for anyone doing similar analysis: a naive lab-frame RMSD (without removing periodic-boundary translation and protein rotation) reports these same stable poses as 8–21 Å "failures"; the protein-frame analysis is what recovers the correct ≈ 2 Å. The full replicated finalist run, including the designed analogs, was in progress when this repository was published and the tables will be updated when it completes.

## Reproduce

The pipeline is packaged as an installable CLI. External tools (docking engine, GNINA, GROMACS, OpenBabel, AmberTools) are system/conda binaries or Docker images and are not installed by pip — see `configs/tools.example.yaml` and `docs/reproducibility_notes.md`.

```bash
conda env create -f environment.yml
conda activate egfr-dockingforge
pip install -e .
pytest -q                       # unit tests
egfrforge --help                # list all stage commands
```

Input structures are public (RCSB PDB; UniProt P00533) and known ligands come from ChEMBL. Each stage records tool versions, command lines, and config hashes in its outputs, and a provenance bundle is produced at the end.

## Repository layout

```text
src/egfr_dockingforge/   pipeline package (stage0 … stage12)
configs/                 one YAML config per stage
workflows/               Snakemake workflow
scripts/                 resumable MD driver and helpers
tests/                   unit tests
docs/                    per-stage protocols, source registries, non-claim notes
```

The Python package keeps its original name `egfr_dockingforge` for import stability; "Syndesis" is the project/method name used in the paper and reports.

## Citation

If you use this work, please cite it via `CITATION.cff`. A manuscript describing the method is in preparation (E. Mitropoulou and D. Giannopoulos, University of Patras). Released under the MIT License.
