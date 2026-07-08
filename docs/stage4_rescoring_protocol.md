# Stage 4 Rescoring Protocol

Stage 4 rescoring evaluates Stage 3 poses without redocking or coordinate refinement. Each pose is scored in place, preserving the Stage 3 receptor, ligand pose, RMSD, sanity, and provisional native-like labels.

The primary scorer is GNINA CNN scoring. GNINA outputs are stored as features, not as proof of binding or activity. Empirical scores from Stage 3 are retained as a classical baseline, and optional Vina/Vinardo rescoring can be added later without changing the table schema.

Required workflow:

1. Build `rescoring_task_matrix.parquet` from Stage 3 poses.
2. Run GNINA in score-only mode for each valid pose.
3. Parse GNINA empirical affinity, CNNscore, CNNaffinity, and CNN_VS when available.
4. Build combined pose score tables and Stage 6 feature exports.
5. Diagnose ranking behavior against Stage 3 RMSD/sanity labels.
6. Generate `reports/04_ml_rescoring.html`.

Stage 4 does not compute ProLIF interactions, train a pose-confidence model, run MD, generate analogs, or screen final candidates.
