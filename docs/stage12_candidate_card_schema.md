# Stage 12 Candidate Card Schema

Each candidate card is a JSON object with stable top-level fields:

- card_version
- final_candidate_id
- molecule_id
- source, subsource, screening_role
- standard_smiles, inchi_key, scaffold_id, novelty_bucket
- closest_known_egfr_ligand
- parent_analog_lineage
- best_pose
- scores
- interactions
- md
- medchem
- evidence_summary
- non_claims
- provenance

Cards are machine-readable evidence records. They are not experimental activity records.
