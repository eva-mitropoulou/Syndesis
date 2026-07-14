# Stage 0 Scope Decisions

Stage 0 fixes the scientific boundary before any docking, rescoring, analog generation,
or MD work is implemented.

## Main Question

Can an interaction-constrained computational workflow optimize noncovalent ATP-site EGFR
analogs while preserving validated EGFR binding modes better than docking-score-only
optimization?

## Included in v1

| Decision | Rationale | Supporting sources |
| --- | --- | --- |
| Human EGFR / ERBB1 / HER1 kinase domain | Keeps the workflow tied to a specific structural target and avoids broad EGFR biology claims. | `klifs_site`, `klifs_2021`, `rcsb_search_api` |
| Experimental holo structures from RCSB PDB and/or KLIFS | Holo structures provide the protein-ligand geometries needed for pose recovery and binding-mode constraints. | `rcsb_search_api`, `klifs_site`, `klifs_2021` |
| Reversible, noncovalent ATP-site small-molecule inhibitors | Keeps the benchmark aligned with pose preservation rather than covalent reaction modeling. | `pdb_1m17`, `pdb_1xkk`, `pdb_4zau` |
| Type I and Type I.5 ATP-site binding modes where classifiable | These binding modes are compatible with ATP-site noncovalent analog optimization. | `klifs_site`, `kincore_site` |
| Active-like and inactive-like receptor states | EGFR has experimentally observed active-like and inactive-like inhibitor-bound structures. | `pdb_1m17`, `pdb_1xkk`, `kincore_site` |
| WT and mutant structures as metadata-bearing records | Mutant structures expand structural evidence, but v1 does not infer mutant selectivity. | `pdb_2ity`, `pdb_2itz` |

## Excluded from v1

| Decision | Rationale | Supporting sources |
| --- | --- | --- |
| Covalent, reversible covalent, and irreversible covalent inhibitors | Covalent inhibition introduces warhead chemistry and reaction geometry outside the v1 objective. | `pdb_4zau` |
| Acrylamide/propiolamide/Michael-acceptor Cys797-targeting warheads | These are covalent-design features rather than noncovalent ATP-site analog features. | `pdb_4zau` |
| Allosteric inhibitors | Allosteric binding is a different pocket and design problem from ATP-site pose preservation. | `klifs_site`, `pdb_1m17`, `pdb_1xkk` |
| Peptides, protein binders, and macrocycles | v1 focuses on conventional small-molecule ATP-site inhibitors. | `rcsb_search_api`, `klifs_site` |
| Fragment-only ligands except optional pocket-mapping references | Fragment poses are useful for mapping but not the initial analog-optimization benchmark. | `posebusters_2023` |
| ATP/ADP/cofactor-only structures except optional references | Cofactors do not represent inhibitor analog binding modes. | `rcsb_search_api` |
| Structures with severe missing active-site residues or unusable ligand identity/geometry | The workflow depends on reliable receptor and ligand coordinates. | `posebusters_2023`, `rcsb_search_api` |

## Validation Boundary

Stage 0 validates project scope, decision logic, and source traceability only. It does not
run docking, GNINA, ProLIF, MD, or deterministic analog generation.
