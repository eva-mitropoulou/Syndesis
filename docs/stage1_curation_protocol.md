# Stage 1 Curation Protocol

Stage 1 constructs the EGFR co-crystal benchmark used by later redocking,
cross-docking, interaction recovery, pose-confidence modeling, and MD validation stages.
It does not implement any of those downstream stages.

## Scope

The benchmark inherits Stage 0 scope:

- Human EGFR / ERBB1 / HER1 kinase domain.
- UniProt accession P00533.
- Orthosteric ATP-binding site.
- Reversible, noncovalent, ATP-site small-molecule inhibitors.
- Active-like and inactive-like receptor states where available.
- WT and mutant structures are allowed, but mutation status is metadata only.

## Acquisition

Candidate PDB IDs are discovered through the RCSB Search API or provided as explicit
control/starter IDs. For each candidate entry, Stage 1 stores:

- mmCIF coordinate file.
- RCSB metadata JSON.
- RCSB validation files when available.
- Chemical Component Dictionary files for native ligands.
- KLIFS metadata when available.

Every downloaded file is tracked with source URL, timestamp, and checksum in acquisition
manifests.

## Complex Definition

One benchmark row represents one PDB entry, one EGFR chain, and one native ligand
instance. The stable complex identifier is:

```text
{pdb_id}_{auth_chain_id}_{ligand_comp_id}_{ligand_instance_id}
```

The workflow does not assume that one PDB entry equals one usable complex.

## Inclusion Rules

Complexes are included when:

- The protein chain maps to human EGFR / ERBB1 kinase domain where metadata are available.
- The structure is experimental.
- The ligand is a small molecule.
- The ligand is bound in or near the ATP pocket.
- The ligand is noncovalent.
- The ligand coordinates are usable.
- ATP/ADP/cofactor-only entries are not treated as primary benchmark complexes.
- Active-site residues are sufficiently modeled.

## Exclusion Rules

Complexes are retained in rejected tables, not hard-deleted, when they violate v1 scope:

- Non-human EGFR.
- Non-kinase-domain receptor.
- No bound small-molecule ligand.
- Ligand outside the ATP pocket.
- Covalent or likely covalent ligand.
- Reversible covalent or irreversible covalent ligand.
- Allosteric ligand.
- Peptide/protein binder.
- Macrocycle in v1.
- Fragment-only ligand in v1 unless reference-only.
- ATP/ADP/cofactor-only ligand unless reference-only.
- Severe missing active-site residues.
- Unusable ligand identity, coordinates, or geometry.

## Native Coordinates

Native ligand and receptor coordinates are exported exactly as deposited. Stage 1 does
not protonate, minimize, charge, or otherwise docking-prepare structures. Later stages
can derive prepared files from these immutable reference exports.

