# Stage 1 Source Registry

| source_id | supported_claim | title | DOI_or_URL | source_type | why_it_matters_for_stage1 |
| --- | --- | --- | --- | --- | --- |
| `rcsb_search_api` | RCSB Search API should be used to search PDB metadata and return PDB identifiers. | RCSB PDB Search API | https://search.rcsb.org/ | documentation | Programmatic candidate structure discovery. |
| `rcsb_data_api` | Detailed RCSB metadata should be fetched through the Data API / GraphQL rather than scraped manually. | RCSB PDB Data API | https://data.rcsb.org/ | documentation | Entry metadata, citation, experimental quality, and entity metadata. |
| `rcsb_mmcif_downloads` | RCSB coordinate downloads provide experimental mmCIF files for parsing native complexes. | RCSB PDB file downloads | https://files.rcsb.org/download/ | data service | Source coordinate files for protein, ligand, and pocket-water extraction. |
| `rcsb_validation_reports` | RCSB validation reports are the preferred source for structure-validation availability. | RCSB validation reports | https://files.rcsb.org/pub/pdb/validation_reports/ | data service | Tracks whether validation XML/PDF files are available. |
| `klifs_database` | KLIFS is appropriate for kinase-specific structure, ligand, and interaction annotation. | KLIFS | https://klifs.net/ | database | Kinase pocket and ligand annotation. |
| `klifs_2021_nar` | KLIFS systematically collects, annotates, and processes kinase structures. | KLIFS overhaul paper | 10.1093/nar/gkaa895 | literature | Source justification for KLIFS metadata. |
| `kincore_database` | Kinase conformational states should be tracked using DFG/C-helix/salt-bridge/HRD/activation-loop metadata. | KinCore | https://dunbrack.fccc.edu/kincore/ | database | Receptor-state metadata policy. |
| `kincore_2019_pnas` | Kinase conformation nomenclature supports structured active/inactive labels. | Defining a new nomenclature for active and inactive kinases | 10.1073/pnas.1814279116 | literature | State classification vocabulary. |
| `egfr_1m17_erlotinib` | EGFR has active-like erlotinib-bound reference structures. | RCSB 1M17 | 10.1074/jbc.M207135200 | structure/literature | Positive active-like reference control. |
| `egfr_1xkk_lapatinib` | EGFR has inactive-like lapatinib-bound reference structures. | RCSB 1XKK | 10.1158/0008-5472.CAN-04-1168 | structure/literature | Positive inactive-like reference control. |
| `egfr_4zau_azd9291_exclusion` | Covalent EGFR inhibitors should be excluded from v1 and used only as negative/scope-exclusion controls. | RCSB 4ZAU | 10.1016/j.jsb.2015.10.018 | structure/literature | Covalent exclusion smoke test. |
| `posebusters_2023` | Native-pose benchmarks should not rely only on RMSD; physical plausibility should be checked later. | PoseBusters | https://arxiv.org/abs/2308.05777 | preprint | Later pose-quality constraints. |
| `interaction_recovery_2024` | Interaction recovery should be part of pose-quality evaluation. | Assessing interaction recovery of predicted protein-ligand poses | https://arxiv.org/abs/2409.20227 | preprint | Later ProLIF-style evaluation rationale. |
| `leakproof_pdbbind_2023` | Similarity/leakage tracking should be planned from benchmark construction. | Leak Proof PDBBind | https://arxiv.org/abs/2308.09639 | preprint | Future split and leakage controls. |

