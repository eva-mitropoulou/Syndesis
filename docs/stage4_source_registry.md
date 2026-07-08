# Stage 4 Source Registry

| source_id | supported_claim | title | DOI_or_URL | source_type | why_it_matters_for_stage4 |
|---|---|---|---|---|---|
| gnina_1_0_2021 | GNINA is an open-source CNN docking/scoring tool. | GNINA 1.0: molecular docking with deep learning | DOI: 10.1186/s13321-021-00522-2 | paper | Supports GNINA as the primary open-source ML rescoring layer. |
| gnina_1_3_2025 | GNINA 1.3 updates the implementation and models. | Gnina 1.3: modern molecular docking with deep learning | https://github.com/gnina/gnina | software | Supports version tracking and modern GNINA usage. |
| ragoza_cnn_scoring_2017 | CNN scoring can score protein-ligand poses. | Protein-Ligand Scoring with Convolutional Neural Networks | DOI: 10.1021/acs.jcim.6b00740 | paper | Foundational CNN scoring source. |
| autodock_vina_2010 | Vina is a classical empirical docking/scoring baseline. | AutoDock Vina: improving the speed and accuracy of docking | DOI: 10.1002/jcc.21334 | paper | Supports original/empirical score comparison. |
| autodock_vina_1_2_2021 | Vina 1.2 supports updated scoring/docking workflows. | AutoDock Vina 1.2.0: new docking methods, expanded force field, and python bindings | DOI: 10.1021/acs.jcim.1c00203 | paper | Supports optional Vina rescoring baseline. |
| posebusters_2024 | Pose evaluation should include physical sanity checks beyond RMSD. | PoseBusters: AI-based docking methods fail to generate physically valid poses | DOI: 10.1039/D3SC04185A | paper | Supports retaining Stage 3 sanity fields in rescoring. |
| interaction_recovery_2024 | Interaction recovery should be evaluated after score-only rescoring. | Assessing interaction recovery of predicted protein-ligand poses | https://arxiv.org/abs/2409.20227 | preprint | Supports deferring interaction-aware labels to Stage 5. |
| fair_comparison_docking_2024 | Docking/scoring comparisons require controlled benchmarking. | Fair comparison of docking methods | https://arxiv.org/ | preprint | Supports scorer-agnostic diagnostics. |
| litpcba_rescoring_2026 | No single scorer dominates all benchmark settings. | LIT-PCBA rescoring benchmark | https://arxiv.org/ | preprint | Supports treating GNINA as feature evidence, not truth. |
| unimol_docking_v2_2024_optional | Modern AI docking can be compared later as a pose-generation method. | Uni-Mol Docking V2 | https://arxiv.org/abs/2405.11769 | preprint | Optional future comparator, not Stage 4 default. |
| diffdock_2022_optional | DiffDock is a generative docking comparator. | DiffDock: Diffusion Steps, Twists, and Turns for Molecular Docking | https://arxiv.org/abs/2210.01776 | preprint | Optional future comparator. |
| nmdn_optional_if_used | Graph scorers may be optional rescoring comparators. | Neural message passing docking/scoring references | https://arxiv.org/ | preprint | Placeholder only if used later. |
