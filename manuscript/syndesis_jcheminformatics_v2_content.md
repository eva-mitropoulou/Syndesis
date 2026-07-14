# Abstract

**Background.** Docking and neural rescoring can rank ligands without establishing whether the selected receptor-specific pose recovers interactions observed in target holo structures. We asked whether a target-native interaction prior adds incremental early-enrichment information to GNINA when both signals are evaluated on the same pose.

**Methods.** The workflow combined four-state EGFR ensemble docking, GNINA CNNscore, and ProLIF interaction fingerprints. A native-derived ATP-site union prior was constructed from the same four EGFR holo complexes. The primary ligand score was the maximum, across those receptor states, of CNNscore multiplied by one plus same-pose interaction recall. The scoring formula and statistical protocol were applied unchanged to a target-specific CDK2 ensemble and native prior. We used paired class-stratified bootstraps, a pose-decoupled late-fusion comparator, three permutation controls, receptor and native-complex exclusions, and overlap, similarity, and size audits.

**Results.** Pose-coupled weighting increased four-receptor EGFR EF1% from 11.98 to 16.40, recovering 89 rather than 65 actives among the first 356 ranked molecules. The paired EF1% difference was 4.42 (95% CI 2.58–6.63). The effect exceeded unrestricted, heavy-atom-count-matched, and class-conditional assignment nulls and persisted after receptor and native-complex exclusions. A pose-decoupled late-fusion comparator reached EF1% 15.85; its difference from pose coupling was unresolved. The target-specific four-receptor/four-native-complex CDK2 analysis was also unresolved (EF1% difference 2.109; 95% CI −0.633 to 4.641). As a downstream stress test, replicated 20 ns MD separated the deliberately mis-docked control from four majority-stable modelled complexes.

**Conclusions.** Native-interaction information improved early enrichment in the evaluated retrospective EGFR benchmark, while same-pose coupling retained an auditable receptor-specific structural basis for the combined score. The unresolved CDK2 result indicates that its contribution must be established separately for each target and receptor ensemble.

**Scientific Contribution.** We introduce and evaluate a pose-coupled native-interaction weighting rule that exposes the structural evidence behind a neural docking score. Paired, permutation, and exclusion controls show that the EGFR early-enrichment signal depends on observed ligand-specific interaction profiles and is not attributable solely to molecular size, activity-class distributions, any single receptor state, or exact native-ligand overlap.

**Keywords:** structure-based virtual screening; docking; interaction fingerprints; early enrichment; EGFR; CDK2

# Background

Structure-based virtual screening is a ranking problem conditioned on pose generation. A docking program may sample a near-native geometry without ranking it first, and a favourable score alone does not establish that the selected pose is structurally credible [@trott2010; @amaro2018; @buttenschoen2024]. Neural scoring functions such as GNINA improve pose assessment by learning protein-ligand patterns from three-dimensional data, but their outputs remain statistical scores rather than explicit tests of whether a pose retains contacts characteristic of a target [@mcnutt2021; @mcnutt2025].

Protein-ligand interaction fingerprints provide that explicit representation. SIFt, SPLIF, kinase interaction profiles, and tools such as ProLIF encode complexes as residue-by-interaction descriptors and have been used for pose comparison, binding-mode analysis, and rescoring [@deng2004; @chuaqui2005; @marcou2007; @da2014; @bouysset2021]. Prior fingerprint-based methods commonly compare a docked pose with a reference pattern after docking. The unresolved methodological issue is how to combine this reference evidence with a learned score across a receptor ensemble.

That combination must preserve pose identity. If the highest neural score is taken from one receptor state and the highest interaction value from another, their combination can reward two incompatible poses. The resulting late-fusion score has no corresponding physical protein-ligand geometry. We instead evaluate both terms for the same ligand, receptor state, and docked pose before maximizing over receptor states.

Syndesis was designed around one narrow hypothesis: a native-derived interaction prior can add early-enrichment information to GNINA when both terms are evaluated for the same receptor-specific pose. EGFR served as the method-development and primary retrospective-evaluation target, whereas CDK2 was used as a target-transfer evaluation after the scoring rule and primary analysis choices had been fixed. The primary endpoint was retrospective ranking rather than activity prediction.

# Methods

## Study design and structural inputs

The primary study was a paired comparison between GNINA and a pose-coupled GNINA-plus-interaction score on the same EGFR ligand-receptor evaluations. The primary EGFR docking ensemble comprised 1M17, 1XKK, 4HJO, and 5CAV. These four ATP-site holo structures also defined the native-derived interaction prior (1M17/AQ4, 1XKK/FMM, 4HJO/AQ4, and 5CAV/4ZQ). A separate five-state sensitivity added ligand-stripped 6DUK while retaining the four-complex ATP-site prior [@to2019]. The primary CDK2 ensemble and native prior both comprised 1FIN/ATP, 2A4L/RRC, 1AQ1/STU, and 1PXN/CK6. Because the extracted 1QMZ receptor lacked deposited phosphothreonine TPO160, it was excluded from the primary CDK2 docking ensemble and native prior. The five-receptor/five-native-complex analysis including this altered representation is reported only as a sensitivity analysis. Receptor identities, chains, docking boxes, quality decisions, and residue maps are machine-readable in the release.

EGFR receptor chains were prepared from PDB structures [@berman2000], with non-protein residues removed and Open Babel 3.1.0 used to produce pH 7.4, Gasteiger-charged PDBQT receptors. Alternate locations were resolved by retaining the selected deposited chain conformer; missing residues and atoms were not modelled. Histidine states followed the Open Babel pH 7.4 protonation assignment. For docked-pose fingerprints, the ProLIF protein was regenerated from the exact docking PDBQT by Open Babel; ProLIF then assigned residue, aromatic, donor, and acceptor chemical perception to that docking-derived model. The EGFR receptor-consistency audit confirmed identical docking and ProLIF heavy-atom sets and zero ProLIF-only atoms for all four primary receptors. CDK2 used the same conversion and fail-closed reconstruction procedure; its receptor preparation and fingerprint status are reported in the release.

**Table 1.** Receptor ensembles, benchmark sizes, and native-prior complexes.

| Target | Docking receptor states | Benchmark molecules | Native-prior complexes |
|---|---|---:|---:|
| EGFR | 1M17, 1XKK, 4HJO, 5CAV | 35,552 (542 actives; 35,010 decoys) | 4 ATP-site complexes |
| CDK2 | 1FIN, 2A4L, 1AQ1, 1PXN | 28,296 (474 actives; 27,822 decoys) | 4 complexes (1FIN/ATP, 2A4L/RRC, 1AQ1/STU, 1PXN/CK6) |

## Docking, interaction encoding, and score coupling

DUD-E input SMILES were converted to one ETKDGv3 three-dimensional state per record (RDKit random seed `0xF00D`), optimised with MMFF94 for at most 1,000 iterations, and converted to PDBQT with Open Babel [@mysinger2012; @landrum2013; @halgren1996; @oboyle2011]. Embedding failures were recorded as preparation failures; an MMFF exception retained the ETKDG conformer and was recorded in the preparation status. Alternative protomers, tautomers, stereoisomers, and conformers were not enumerated. Each ligand was docked independently to the target-specific receptor ensemble with Uni-Dock 1.2.0 in `balance` mode, nine output modes, and seed 807 [@yu2023unidock]. The four CDK2 receptors used the same Uni-Dock, GNINA, ProLIF, molecular-state, and graph-reconstruction settings as EGFR. Only the top Uni-Dock pose per ligand-receptor pair entered the enrichment analysis. GNINA 1.3.3 used its default `rescore` CNN model in `--score_only` mode, with pose minimisation disabled; the parsed `CNNscore` field was the neural ranking term [@mcnutt2021; @mcnutt2025].

ProLIF 2.2.0 encoded hydrophobic, implicit donor and acceptor, ionic, cation-pi, pi-cation, and van der Waals interactions as normalised residue-by-interaction-type bits [@bouysset2021]. Configured distance cutoffs were 3.6 Å for hydrogen bonds, 4.5 Å for hydrophobic and ionic contacts, and 4.0 Å for van der Waals contacts; ProLIF default geometry policies otherwise applied. Docked coordinates were transferred onto the prepared SDF graph before fingerprinting, preserving bond order, formal charge, tautomerism, and stereochemistry. The mapping matched each prepared-SDF heavy atom to the pre-docking PDBQT heavy atom by element and coordinate within 0.05 Å, required exact prepared/posed atom count and element order, then transferred the corresponding docked coordinate to the unchanged RDKit graph. Graph identity, formal charges, and stereocentres therefore came from the prepared SDF rather than PDBQT. Atom-order, element, count, and coordinate-mapping failures were not converted to zero scores; all 142,208 primary EGFR ligand-receptor evaluations passed this reconstruction audit.

Let $F_{i,r}$ be the fingerprint for ligand $i$ in receptor state $r$, and $N_k$ the fingerprint for native complex $k$. The primary prior was the target-native union $C=\bigcup_k N_k$, with 62 EGFR and 38 CDK2 bits. Same-pose recall was $R_{i,r}=|F_{i,r}\cap C|/|C|$. The coupled ligand score was

$$
S_i=\max_r\{\mathrm{CNNscore}_{i,r}[1+R_{i,r}]\}.
$$

The matched GNINA baseline was $G_i=\max_r\{\mathrm{CNNscore}_{i,r}\}$ over the same four receptor-specific top Uni-Dock poses. The pose-decoupled comparator was $L_i=[\max_r\{\mathrm{CNNscore}_{i,r}\}][1+\max_r\{R_{i,r}\}]$; it can combine values from different receptor states and is therefore not assigned to a single physical pose.

The multiplier is bounded between one and two. Critically, CNNscore and recall in each product belong to the same pose. As summarized in @fig-pose-coupling-workflow, the workflow preserves this same-pose relationship. EGFR labels informed method development, so EGFR is reported as the method-development and primary retrospective-evaluation target. The scoring formula (multiplicative union recall with $\lambda=1$), EF1% endpoint, bootstrap design, and seeds were fixed before final strict fingerprint recomputation and applied unchanged to the target-specific CDK2 receptor ensemble and native prior. Alternative priors and formulas were sensitivity analyses, not selection criteria.

![Same-pose interaction-aware ranking workflow. Four ATP-site EGFR native complexes define a union interaction prior. Each candidate is docked independently into the four primary receptor states; CNNscore and recall are evaluated for the same selected pose before the maximum coupled score is used for ligand ranking.](figures/figure1_pose_coupling_workflow.png){#fig-pose-coupling-workflow width=72% fig-alt="Vertical workflow from four native EGFR complexes to same-pose coupled ligand ranking across four receptor states."}

## Statistical evaluation

EF1% was the primary endpoint; ROC-AUC, EF5%, and BEDROC ($\alpha=80.5$) were secondary [@truchon2007]. The top 1% set used $\max[1,\operatorname{round}(0.01N)]$ molecules, with stable input order resolving score ties. We used 2,000 paired class-stratified bootstrap resamples (seed 807). Three 1,000-draw permutation controls reassigned complete receptor-ensemble recall vectors: across all ligands, within heavy-atom-count deciles, and within activity class. The last preserves active-decoy recall distributions and tests molecule-specific assignment rather than a general random-score null. Empirical permutation values were $p=(b+1)/(1000+1)$, where $b$ was the number of null draws at least as large as the observed EF1%. Receptor exclusions, native-complex exclusions, joint duplicate-chemotype exclusions, exact native/DUD-E identity checks, ECFP4 similarity strata, and recall-size correlations tested robustness [@bemis1996; @rogers2010]. ECFP4 used RDKit Morgan fingerprints with radius 2 and 2,048 bits, without chirality or feature invariants. The five-receptor 6DUK-inclusive result was evaluated separately as an ensemble sensitivity, and receptor-specific-prior sensitivity was restricted to the four primary EGFR receptors. DUD-E decoys are property-matched benchmark compounds rather than experimentally confirmed inactives, and its analogue and decoy construction can introduce benchmark-specific bias; this analysis therefore evaluates incremental ranking performance, not prospective activity prediction [@mysinger2012; @stein2021; @wallach2018].

## Redocking evaluation

The manuscript redocking check used the four ATP-site cognate tasks 1M17/AQ4, 1XKK/FMM, 4HJO/AQ4, and 5CAV/4ZQ. Uni-Dock 1.2.0 ran in `balance` mode with exhaustiveness 8, 20 output modes, and seed 13. Rectangular boxes were derived from each cognate deposited ligand; the released redocking-task manifest records each box centre and side length. Native-like poses were defined by in-place, symmetry-corrected heavy-atom RMSD of at most 2.0 Å to the receptor-frame deposited ligand. RMSD used RDKit template-based bond-order assignment and `CalcRMS` without coordinate superposition; failures were retained as failed mappings rather than silently substituted.

## Molecular-dynamics stress test

Molecular dynamics (MD) was a downstream pose-persistence stress test, not an activity or binding-free-energy calculation. Seven preselected 1M17 EGFR complexes were simulated: three known-ligand controls, three deterministic RDKit-generated analogues, and one deliberately mis-docked negative control. The three controls were included as known EGFR ligands; the three analogues were the highest-priority accepted products of deterministic rule-based transformations of the Control 002 parent, selected before trajectory analysis because they preserved the selected binding mode without a major score or ligand-efficiency loss. The negative control used the Control 002 ligand translated by 8.0 Å and rotated by 180 degrees from its selected pose. The released system-selection manifest records all seven systems; no system was removed after preparation or trajectory inspection.

**Table 2.** MD-system selection and provenance. All systems used receptor 1M17 and the selected pose 01. The released `md_candidate_manifest.csv` provides the exact prepared SMILES, pose identifiers, and input paths for every row.

| System | Ligand identifier | Role | Parent and deterministic rule | Pre-MD selection reason |
|---|---|---|---|---|
| Control 001 | `mol_aatoplajlqecrd` | Known-ligand control | Self | Predeclared known-ligand control |
| Control 002 | `mol_aakjlrggtjkamg` | Known-ligand control | Self | Predeclared known-ligand control |
| Control 003 | `mol_aaalvybiclmama` | Known-ligand control | Self | Predeclared known-ligand control |
| Analogue 004 | `analog_e1bd229bfc1e` | Deterministic analogue | Control 002; small-substituent scan | Tier 2, preserved binding mode |
| Analogue 005 | `analog_0367d4e7eae4` | Deterministic analogue | Control 002; halogen scan | Tier 2, preserved binding mode |
| Analogue 006 | `analog_2939c889a4f4` | Deterministic analogue | Control 002; halogen scan | Tier 2, preserved binding mode |
| Mis-docked control | Control 002 | Negative control | +8.0 Å, 180-degree displacement | Deliberate gate challenge |

Ligands retained their single standardised input molecular state and were parameterised with AmberTools 24.8 using the AmberTools GAFF2 implementation and AM1-BCC charges, then exported through ACPYPE 2023.10.27 [@case2023; @jakalian2000; @jakalian2002; @daSilva2012]. Each complex used Amber ff19SB protein parameters and OPC3 water [@tian2020; @izadi2016; @abraham2015]. Systems used a dodecahedral box with 1.0 nm padding, were neutralised, and were brought to 0.15 M NaCl.

Each system underwent energy minimisation (up to 50,000 steps), 0.5 ns restrained NVT and 1.0 ns restrained NPT equilibration at 300 K and 1 bar, followed by three independent 20 ns NPT production trajectories with a 2 fs timestep. All protein heavy atoms were restrained during equilibration with force constants of 1,000 kJ mol$^{-1}$ nm$^{-2}$ in each Cartesian direction; no ligand restraints were enabled. GROMACS 2026.0 used LINCS constraints on covalent bonds involving hydrogen atoms, a V-rescale thermostat for protein and non-protein groups (ligand, water, and ions; $\tau_T=0.1$ ps), isotropic Parrinello--Rahman pressure coupling ($\tau_P=2.0$ ps, reference pressure 1 bar, compressibility $4.5\times10^{-5}$ bar$^{-1}$), particle-mesh Ewald electrostatics with a 1.2 nm real-space cutoff, and a 1.2 nm Lennard-Jones cutoff. The archived TPR files retain resolved Verlet neighbour-list settings; frames were saved every 10 ps. Replicates 2 and 3 started from new Maxwell velocities at 300 K with seeds 20262708 and 20263708; replicate 1 continued from equilibration. No production burn-in was discarded. MDAnalysis 2.9.0 reconstructed trajectories across periodic boundaries and aligned each frame to the initial backbone atoms of pocket residues within 6 Å of the reference ligand. Minimum-image distances were used for all distance-based occupancies. Ligand heavy-atom RMSD was then measured in that local protein frame. The predeclared *pocket-contact retention* gate required any ligand heavy atom to be within 4.5 Å of the starting-pocket heavy-atom set in at least 90% of frames; ligand centre-of-mass drift was retained as a complementary, stricter displacement measure and is reported with the full metrics. This permissive proximity criterion is not interpreted as a stand-alone binding measure.

The MD interaction set was a preconfigured conserved EGFR core: residue--interaction bits observed in at least 60% of the four native complexes, plus the three overrides fixed in `configs/stage5_interaction_atlas.yaml` before MD (hinge HBAcceptor, hinge HBDonor, and gatekeeper Hydrophobic). Hinge occupancy was the mean occupancy of the core bit(s) assigned the `hinge` role; key-contact occupancy was the mean across all core bits. Because the trajectory implementation evaluates minimum ligand--residue heavy-atom distances, these are reported as *distance-based contact occupancies*, including for bits originally typed as hydrogen-bond interactions; they are not directional hydrogen-bond geometries. Cutoffs were 3.5 Å for hydrogen-bond-typed contacts, 4.0 Å for ionic contacts, and 4.5 Å otherwise. The full per-frame measurements and interaction identities are released.

A replicate was labeled stable only when median ligand RMSD was at most 3.0 Å, its 95th-percentile RMSD was at most 5.0 Å, at least 90% of frames remained in the pocket, hinge occupancy was at least 0.30, and mean key-interaction occupancy was at least 0.50. A system required at least two of its three replicates to be stable. Parameterization and trajectory warnings were retained in the released tables; GAFF2/AM1-BCC parameters are an open, practical parameterization workflow rather than a guarantee of exact ligand physics.

# Results

## Redocking illustrates distinct pose-generation and pose-selection outcomes

Before evaluating ligand enrichment, we performed redocking as a methodological check to determine whether Uni-Dock could generate native-like poses and whether its scoring function selected them. In the two AQ4 tasks, 1M17 and 4HJO sampled strict native-like poses at 1.63 Å (rank 3) and 1.43 Å (rank 4), whereas their Uni-Dock rank-1 poses were 7.95 Å and 6.67 Å from the deposited ligands. The 1XKK/FMM task sampled a 1.52 Å pose at rank 3, while 5CAV/4ZQ selected a 0.32 Å pose at rank 1. This four-task check evaluates pose sampling and docking pose selection only; it does not test the interaction-coupled enrichment rule or estimate a general redocking ranking metric.

## Pose-coupled weighting improves early EGFR enrichment

Table 3 and @fig-enrichment summarize the primary comparison between GNINA and the pose-coupled score on the four-receptor EGFR ensemble. The benchmark contained 542 actives and 35,010 decoys, and the top 1% of the ranking therefore comprised 356 molecules. GNINA recovered 65 actives in this subset, corresponding to an EF1% of 11.98, whereas the pose-coupled score recovered 89 actives and increased EF1% to 16.40. This represents 24 additional actives and a paired EF1% improvement of 4.42 units (95% CI 2.58–6.63). Paired effects also favoured coupling for ROC-AUC (0.00497, 95% CI 0.00362–0.00639), EF5% (0.70, 0.30–1.07), and BEDROC (0.072, 0.054–0.090). The largest practical effect occurred at the top of the ranking, consistent with the intended use of the method for early compound prioritization.

**Table 3.** EGFR enrichment across the four-receptor primary ensemble. Intervals are percentile 95% confidence intervals from 2,000 class-stratified bootstrap resamples.

| Ranking arm | ROC-AUC (95% CI) | EF1% (95% CI) | EF5% (95% CI) | BEDROC (95% CI) |
|---|---:|---:|---:|---:|
| GNINA | 0.770 (0.746–0.794) | 11.98 (9.40–14.56) | 7.01 (6.27–7.78) | 0.210 (0.178–0.244) |
| **Pose-coupled score** | **0.775 (0.751–0.798)** | **16.40 (13.63–19.35)** | **7.71 (6.90–8.52)** | **0.282 (0.245–0.320)** |

The corresponding paired effects (pose-coupled minus GNINA) were +0.00497 ROC-AUC (95% CI 0.00362–0.00639), +0.70 EF5% (0.30–1.07), and +0.0721 BEDROC (0.0544–0.0898).

![EGFR primary result and CDK2 transfer boundary. Error bars are class-stratified bootstrap 95% confidence intervals.](figures/figure2_enrichment.png){#fig-enrichment width=100% fig-alt="Enrichment metrics for GNINA and pose-coupled scoring on EGFR and CDK2."}

We next tested whether the EF1% increase could be reproduced after disconnecting interaction recall from the ligand to which it belonged. Each ligand’s complete four-receptor recall vector was randomly reassigned 1,000 times under three permutation schemes. The unrestricted permutation exchanged recall vectors among all ligands and produced a mean EF1% of 11.35 ($p=0.0010$). The heavy-atom-count-matched permutation exchanged vectors only among similarly sized ligands and produced a mean EF1% of 12.39 ($p=0.0010$), indicating that the observed gain was not explained by ligand size alone. The class-conditional permutation exchanged vectors separately among actives and among decoys, thereby preserving class-level recall differences while disrupting the association between each ligand and its own interaction profile. This more stringent null produced a mean EF1% of 14.26 ($p=0.0040$). As shown in @fig-permutation, the observed EF1% of 16.40 lay beyond all three permutation distributions, supporting a molecule-specific contribution from the correctly assigned interaction information.

**Table 4.** Pose-decoupled late-fusion sensitivity on the primary EGFR ensemble. Point estimates are shown for ranking arms; paired effects are pose-coupled minus pose-decoupled late fusion and use class-stratified bootstrap 95% confidence intervals.

| Contrast or arm | ROC-AUC | EF1% | EF5% | BEDROC |
|---|---:|---:|---:|---:|
| Pose-decoupled late fusion | 0.776 | 15.85 | 7.60 | 0.269 |
| Pose-coupled minus late fusion | −0.00106 (−0.00217 to −0.000061) | 0.55 (−0.55 to 2.03) | 0.11 (−0.11 to 0.44) | 0.014 (0.006 to 0.022) |

The pose-decoupled late-fusion comparator independently combined each ligand’s highest CNNscore and highest recall across receptor states, even when the two values originated from different poses. It achieved an EF1% of 15.85 and recovered 86 actives. The paired comparison did not establish superior EF1% enrichment for pose coupling: late fusion had slightly higher ROC-AUC, whereas pose coupling had higher BEDROC; EF5% superiority was likewise unresolved. Its structural advantage is instead that the neural and interaction components can be attributed to one receptor-specific ligand pose. In 21,993 of 35,552 EGFR ligands (61.9%), the receptor state with maximum CNNscore differed from that with maximum recall, quantifying the potential structural mismatch avoided by pose coupling.

![EF1% permutation distributions for EGFR and CDK2 under unrestricted, heavy-atom-count-matched, and class-conditional interaction assignments.](figures/figure4_permutation_nulls.png){#fig-permutation width=100% fig-alt="Permutation distributions with observed pose-coupled enrichment marked for EGFR and CDK2."}

## Robustness and ensemble sensitivity of the EGFR result

The EGFR enrichment gain was not driven by any single receptor state. When each of the four primary receptors was removed in turn, the paired EF1% improvement remained positive, ranging from 3.50 to 4.79, and all 95% confidence intervals excluded zero (@fig-egfr-receptor-sensitivity). Potential overlap between the native ligands and the benchmark also did not explain the result. Rebuilding the native-union prior without 1XKK/FMM, while retaining the 1XKK docking receptor and all benchmark molecules, retained an EF1% of 15.29 and a paired gain of 3.32 (95% CI 1.84–5.34). Rebuilding the prior without both AQ4-containing native complexes, while retaining their receptor states and all benchmark molecules, retained an EF1% of 15.66 and a gain of 3.69 (1.11–6.08). Moreover, among the 369 actives with maximum ECFP4 similarity below 0.30 to every distinct native ligand, the coupled ranking recovered 54 in the global top 1%, compared with 39 for GNINA.

![Leave-one-receptor-out robustness of EGFR early enrichment. Each cell shows the paired EF1% difference between the pose-coupled score and GNINA after exclusion of one receptor from the four-receptor primary ensemble. The enrichment gain remained positive after every receptor exclusion, indicating that no single receptor state accounted for the primary EGFR result.](figures/figure5_receptor_sensitivity.png){#fig-egfr-receptor-sensitivity width=86% fig-alt="One-row green heatmap of positive paired EF1 percent differences after each primary EGFR receptor exclusion."}

Native-union recall was weakly associated with heavy-atom count (Spearman $\rho=0.218$) and molecular weight ($\rho=0.221$), but more strongly associated with the total number of detected contacts ($\rho=0.704$). Recall should therefore not be interpreted as independent of ligand size or contact abundance. Nevertheless, the heavy-atom-count-matched permutation analysis showed that molecular size alone did not account for the observed enrichment gain.

**Table 5.** Alternative EGFR interaction-prior sensitivities. Effects are paired EF1% differences from GNINA with 2,000 class-stratified bootstrap 95% confidence intervals.

| Prior definition | EF1% | Paired difference from GNINA (95% CI) |
|---|---:|---:|
| 60% conserved core | 17.14 | 5.16 (3.13–7.37) |
| Frequency-weighted recall | 16.58 | 4.61 (2.95–6.82) |
| Jaccard similarity | 16.03 | 4.05 (2.21–6.45) |
| Receptor-specific recall | 16.95 | 4.97 (3.13–7.55) |

The positive result was retained under alternative definitions of the interaction prior. The conserved core retained bits found in at least 60% of native complexes; frequency weighting used each bit's native-complex frequency; Jaccard used $|F\cap C|/|F\cup C|$; and receptor-specific recall used the matching receptor's native fingerprint. Across the multiplicative $\lambda$ sweep, EF1% remained above the GNINA baseline from $\lambda=0.25$ to 3, ranging from 13.82 to 18.24. These analyses support multiple viable formulations, with union recall at $\lambda=1$ retained as the development-fixed, interpretable primary definition.

Ensemble composition was examined separately by adding ligand-stripped 6DUK to the four primary receptor states while retaining the same ATP-site interaction prior. GNINA EF1% changed slightly from 11.98 to 11.79, whereas the coupled EF1% remained 16.40 with the same 89 actives recovered in the top 1%. The five-receptor paired improvement was 4.61 (95% CI 2.58–6.82). Thus, adding 6DUK did not materially alter the coupled ranking or the study conclusion; 6DUK was retained only as an ensemble-sensitivity analysis and was not part of the primary protocol.

## CDK2 defines a transfer boundary

CDK2 was used to test transfer of the EGFR-developed scoring formula to a second kinase target. Across the primary four-receptor/four-native-complex analysis (1FIN, 2A4L, 1AQ1, and 1PXN), coupled scoring increased EF1% from 11.39 to 13.50, corresponding to a paired difference of 2.109 (95% CI −0.633 to 4.641; bootstrap proportion above zero 0.928). The paired interval did not exclude zero, so this positive result remains unresolved. The observed EF1% also exceeded the unrestricted, heavy-atom-count-matched, and class-conditional assignment nulls: means 9.30, 11.70, and 11.73, observed-minus-null differences 4.20, 1.80, and 1.77, and empirical $p=0.0010$, 0.0300, and 0.0490, respectively. These permutation tests assess whether observed ligand--interaction assignments outperform randomised assignments, whereas the paired bootstrap assesses whether coupling reliably outperforms GNINA. The five-receptor/five-native-complex analysis including altered 1QMZ is retained only as a sensitivity analysis. Native-prior overlap also remained unresolved: rebuilding the four-complex prior without the remaining ATP complex 1FIN/ATP gave EF1% 13.08 (paired difference 1.69; 95% CI −1.05 to 4.01), while rebuilding it without exact-overlap inhibitor complexes 2A4L/RRC and 1AQ1/STU gave EF1% 13.08 (difference 1.69; −0.63 to 4.64). These sensitivity results reinforce that CDK2 transfer should be interpreted conservatively.

**Table 6.** CDK2 four-receptor/four-native-complex transfer analysis. Intervals are percentile 95% confidence intervals from 2,000 class-stratified bootstrap resamples.

| Ranking arm | ROC-AUC (95% CI) | EF1% (95% CI) | EF5% (95% CI) | BEDROC (95% CI) |
|---|---:|---:|---:|---:|
| GNINA | 0.749 (0.721–0.776) | 11.39 (9.07–13.92) | 6.88 (6.03–7.68) | 0.229 (0.192–0.264) |
| Pose-coupled score | 0.753 (0.725–0.780) | 13.50 (10.55–16.24) | 7.05 (6.20–7.89) | 0.248 (0.209–0.287) |

![CDK2 receptor-exclusion sensitivity. Points show paired EF1% differences and horizontal bars show paired 95% confidence intervals after exclusion of one receptor from the four-receptor CDK2 transfer ensemble.](figures/figure6_cdk2_receptor_sensitivity.png){#fig-cdk2-receptor-sensitivity width=78% fig-alt="Forest plot of CDK2 receptor-exclusion paired EF1 percent differences and 95 percent confidence intervals."}

## Replicated MD distinguishes persistent and unstable poses

All 21 planned production trajectories completed. The majority-replicate gate labelled two of three known-ligand controls and two of three deterministic analogues as MD-stable; the deliberately mis-docked control was unstable in all three replicates. The latter had a median ligand RMSD of 5.72 Å and median key-contact occupancy of 0.008, whereas accepted systems had median ligand RMSD of 1.80–2.30 Å and key-contact occupancy of 0.52–0.69. Control 001 was rejected in all three replicates despite a median RMSD of 1.99 Å because it retained the permissive pocket criterion but lost the hinge contact (occupancy 0.00–0.02) and mean core-contact occupancy (0.41–0.44); the rejection therefore reflects interaction loss from the docked starting pose, not ligand egress. Analogue 005 failed geometric criteria in two replicates and contact criteria in the remaining replicate. These labels assess persistence of the modelled binding mode and do not establish affinity or biological activity.

**Table 7.** Replicated 20 ns MD gate summary. Values are candidate-level medians except pocket retention, which is the minimum across replicates. The majority decision was computed from complete replicate-level values, not from these summaries alone. Ligand identities and provenance are in Table 2.

| System | Stable reps | RMSD med./p95 (Å) | Pocket-contact | Hinge/key occupancy | Failed gate | Final label |
|---|---:|---:|---:|---:|---|---|
| Control 001 | 0/3 | 1.99 / 3.53 | 1.00 | 0.007 / 0.433 | Contacts, all reps | Unstable |
| Control 002 | 3/3 | 1.80 / 2.44 | 1.00 | 0.774 / 0.618 | None | Stable |
| Control 003 | 3/3 | 2.16 / 2.54 | 1.00 | 0.997 / 0.688 | None | Stable |
| Analogue 004 | 2/3 | 1.87 / 3.11 | 1.00 | 0.905 / 0.519 | Contacts, one rep | Stable |
| Analogue 005 | 0/3 | 3.43 / 4.07 | 1.00 | 0.941 / 0.550 | Geometry, two reps; contacts, one rep | Unstable |
| Analogue 006 | 3/3 | 2.30 / 3.00 | 1.00 | 0.472 / 0.644 | None | Stable |
| Mis-docked control | 0/3 | 5.72 / 8.43 | 1.00 | 0.000 / 0.008 | Geometry, all reps | Unstable |

![MD pose-persistence stress test. Points summarise candidate-level median ligand RMSD with the range across three independent 20 ns replicates; bars show candidate-level median core-contact occupancy. Colours indicate the predeclared majority-replicate decision. The deliberately mis-docked control is separated from the majority-stable systems by both geometric and interaction evidence. Full-protein start/end coordinate renders for the stable and failed examples are retained in the public repository.](figures/figure6_md_stability.png){#fig-md width=100% fig-alt="Candidate-level ligand RMSD and core-contact occupancy for replicated MD systems, including a deliberately mis-docked negative control."}

# Discussion

The main finding is methodological and specific: a native-derived interaction prior improved EGFR early enrichment when it was coupled to the GNINA score of the same receptor-specific pose. Same-pose coupling preserves a receptor-specific structural origin for the combined score because an ensemble otherwise permits components from incompatible poses to be combined. The permutation analyses show that the enrichment gain depends on observed ligand-specific interaction profiles rather than arbitrary, size-matched, or activity-class-preserving assignments. The pose-decoupled late-fusion rule achieved similar EF1%, however, and the paired comparison did not establish a performance advantage for same-pose coupling. Its demonstrated advantage is therefore interpretability: the combined score remains attributable to one receptor-specific protein-ligand geometry.

Previous interaction-fingerprint methods have shown that native patterns can aid pose comparison and ranking. The present design quantifies the incremental contribution of target-native interaction information while preserving the receptor-specific structural origin of the combined score. The robustness analyses refine the biological interpretation: the effect survived removal of individual receptor states, exact native-ligand overlap, and duplicate AQ4 structures, and remained visible among low-similarity actives. These results support the view that the prior contributes target-structural information rather than merely recognizing one crystallographic chemotype. At the same time, union recall depends on the number of detected contacts and different prior definitions also perform well. The sensitivity analyses indicate that several formulations of target-native interaction evidence can complement GNINA, while union recall provides a simple and interpretable primary definition.

CDK2 sets the boundary of the enrichment claim. Its point estimates are favourable, but the paired EF1% interval did not exclude zero and receptor dependence is substantial. Transfer should therefore be evaluated target by target rather than assumed from kinase-family membership. The MD stress test adds a different type of evidence: it can distinguish persistent modelled poses from a deliberately mis-docked control, but it does not validate retrospective enrichment, calculate binding free energies, or establish experimental activity. Because only the top Uni-Dock mode for each ligand-receptor pair was rescored and fingerprinted, the present analysis evaluates ligand ranking and receptor-state selection; it does not test whether interaction weighting can recover lower-ranked modes sampled within a receptor. The one-state ligand-preparation strategy is an additional practical limitation because alternative protonation, tautomeric, stereochemical, and conformational states can affect both docking and interaction recovery. Likewise, DUD-E is external to Syndesis development but not necessarily independent of GNINA training structures or chemistry; the paired design estimates the incremental contribution of interaction coupling within the evaluated benchmarks, while broader generalisation requires benchmarks independent of both method development and neural-score training. The primary EGFR ensemble was restricted to ATP-site holo conformations; adding the allosteric-ligand-stabilised 6DUK state was examined separately to avoid conflating primary receptor selection with a structurally distinct sensitivity.

# Conclusions

Pose-coupled native-interaction weighting improved EGFR early enrichment beyond GNINA across paired, permutation, four-primary-receptor exclusion, native-overlap, and similarity controls. A pose-decoupled fusion rule achieved comparable enrichment, indicating that the specific empirical ranking advantage of same-pose coupling remains unresolved; its methodological advantage is that the combined score remains attributable to one receptor-specific protein-ligand geometry. The unresolved CDK2 transfer result further indicates that native-interaction priors should be constructed and validated at the target-and-ensemble level.

# Abbreviations

ATP, adenosine triphosphate; BEDROC, Boltzmann-enhanced discrimination of receiver operating characteristic; CDK2, cyclin-dependent kinase 2; CNN, convolutional neural network; DUD-E, Directory of Useful Decoys--Enhanced; ECFP, extended-connectivity fingerprint; EF, enrichment factor; EGFR, epidermal growth factor receptor; IFP, interaction fingerprint; MD, molecular dynamics; NDCG, normalised discounted cumulative gain; ROC-AUC, area under the receiver-operating-characteristic curve.

# Declarations

## Availability of data and materials

Source code, workflow configurations, tests, figures, the rendered manuscript, and machine-readable supporting data are available at [https://github.com/eva-mitropoulou/Syndesis](https://github.com/eva-mitropoulou/Syndesis). The exact paper package is the [`v1.1.7-paper`](https://github.com/eva-mitropoulou/Syndesis/tree/v1.1.7-paper) tagged repository snapshot. It includes pose coordinates, native interaction-bit tables, ligand-level benchmark scores and fingerprints, bootstrap and permutation draws, late-fusion and exclusion analyses, graph-mapping validation, four-receptor primary manifests, AMBER/GAFF2 parameterization reports, GROMACS inputs, replicate-level MD metrics, and per-frame pose-persistence measurements. Raw structures and benchmark molecules originate from the PDB and DUD-E and remain subject to their source terms. No separate supplementary document accompanies this manuscript.

## Competing interests

The authors declare no competing interests.

## Generative-AI assistance

Generative-AI tools were used for language editing, code-drafting and review assistance, and literature scoping. The authors reviewed and validated all code, analyses, figures, citations, and claims, and retain full responsibility for the manuscript.

## Funding

This research received no specific grant from any funding agency in the public, commercial, or not-for-profit sectors.

## Authors' contributions

Following the CRediT taxonomy: E.M., conceptualization, methodology, software, formal analysis, investigation, data curation, visualization, writing--original draft, and writing--review and editing; D.G., conceptualization, methodology, validation, resources, supervision, project administration, and writing--review and editing. Both authors read and approved the manuscript.

## Ethics approval and consent to participate

Not applicable.

## Consent for publication

Not applicable.

## Acknowledgements

The authors acknowledge the developers and maintainers of the open scientific software and public structural and chemical databases used in this study.
