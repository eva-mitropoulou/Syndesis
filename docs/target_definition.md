# Target Definition

## Target

Syndesis is a kinase-docking methodology evaluated here on the human epidermal growth
factor receptor kinase domain as its primary case study. CDK2 provides a separate
target-transfer analysis. The EGFR case-study scope is:

- Primary target: human EGFR / ERBB1 / HER1 kinase domain
- Organism: Homo sapiens
- UniProt accession: P00533
- Binding site: orthosteric ATP site
- Ligand class in v1: reversible, noncovalent, ATP-site small-molecule inhibitors

The project is a computational-only structure-based design workflow. It is not a broad
EGFR drug-discovery program.

## Why EGFR Kinase Domain Instead of Broad EGFR Drug Discovery

EGFR is a multidomain receptor tyrosine kinase, but this project is limited to the
intracellular catalytic kinase domain because the workflow depends on experimentally
resolved protein-ligand binding geometries. The project question is about preserving
validated ATP-site binding modes during analog optimization, so extracellular-domain
biology, antibody binding, receptor dimerization, signaling assays, and cellular potency
are outside v1.

This narrower definition makes the benchmark testable: the workflow can compare docked
poses against kinase-domain holo structures, residue-level interaction annotations, and
kinase-state metadata from structure-focused resources such as RCSB PDB, KLIFS, and
KinCore.

## Why v1 Is Noncovalent ATP-Site Only

EGFR has reversible ATP-site inhibitors and covalent Cys797-targeting inhibitors. These
are different modeling problems. Covalent design requires reaction geometry, warhead
reactivity, covalent bond formation, residue protonation/tautomer assumptions, and
mechanism-specific validation. Including those compounds in v1 would blur the project
claim and make docking-score comparisons less interpretable.

Therefore v1 includes reversible, noncovalent ATP-site ligands only. Type I and Type I.5
ATP-site binding modes are included where they can be classified. Covalent inhibitors,
reversible covalent inhibitors, irreversible covalent inhibitors, and Cys797-targeting
Michael-acceptor warheads are excluded.

## Why Experimental Holo Structures Are Primary

The core evidence for this project is known EGFR binding geometry. Experimental holo
structures provide ligand coordinates, protein conformations, and residue-level contacts
that can be used to define what credible EGFR binding looks like. RCSB PDB provides
experimental coordinate entries, while KLIFS adds kinase-focused pocket and interaction
annotations.

Computed apo structures are not the primary source for v1 because the workflow needs
ligand-bound receptor conformations. AlphaFold and related computed models can be useful
fallbacks for missing regions or exploratory comparisons, but they are not treated as
substitutes for experimental EGFR holo structures when those structures are available.

## Mutation Policy

WT and mutant EGFR structures are both allowed in the structural benchmark because both
can inform receptor-state diversity and pose recovery. However, mutation status is tracked
as metadata only in v1.

The project will not claim mutant selectivity, resistance coverage, or genotype-specific
activity. Any future mutant-selectivity claim would require a separate design objective,
balanced WT/mutant benchmarks, activity data, and validation metrics.

## Receptor-State Policy

EGFR ATP-site inhibitors bind more than one kinase conformation. Active-like and
inactive-like receptor states should both be included when supported by experimental holo
structures. The project tracks structural metadata such as DFG state, C-helix/salt-bridge
state, activation-loop state, mutation status, ligand identity, and resolution rather than
depending on vague active/inactive labels alone.

This state diversity is necessary because a workflow trained on one receptor state may
overfit to one binding geometry and reject plausible alternatives.

## Claim Boundary

The project may claim:

- a reproducible computational workflow for structure-constrained EGFR analog optimization
- a curated EGFR kinase-domain holo-structure benchmark
- pose and analog selection based on docking scores plus interaction and pose-confidence constraints
- comparison against docking-score-only optimization

The project will not claim:

- wet-lab validation
- experimentally confirmed inhibitor discovery
- mutant-selective inhibitor design in v1
- covalent inhibitor design in v1
- allosteric inhibitor design in v1
- clinical, cellular, or biochemical activity unless external data are explicitly modeled later
