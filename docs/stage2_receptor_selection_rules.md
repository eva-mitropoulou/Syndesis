# Stage 2 Receptor Selection Rules

Candidate receptors must pass Stage 1 inclusion and Stage 2 hard filters:

- `include_in_stage1_benchmark == true`
- `covalent_flag == false`
- `allosteric_flag == false`
- `atp_site_flag == true`
- receptor and native ligand files exist
- active-site completeness score passes the configured minimum
- quality tier is allowed by config
- resolution passes the configured maximum

Within state strata, receptors are clustered by pocket geometry. Medoids are selected by
cluster centrality and quality. Reference controls are considered if they pass filters,
but controls are not blindly forced into the ensemble if they fail scope or quality.

