# Stage 4 Score Interpretation

Stage 4 compares score ranking directions carefully:

- original docking score: lower is better
- Vina/Vinardo empirical scores: lower is better
- GNINA empirical affinity: lower is better
- CNNscore: higher is better
- CNNaffinity: higher is better for ranking in this workflow
- CNN_VS: higher is better when available or derived

Native-like labels are inherited from Stage 3 and remain provisional. Interaction recovery is deliberately marked `pending_stage5`; final pose labels are also `pending_stage5`.

Useful Stage 4 outcomes include:

- GNINA rescues a docking-score ranking failure
- GNINA prefers a non-native-like pose
- GNINA prefers an invalid pose
- empirical and CNN scores strongly disagree
- scoring fails for a receptor/ligand state

These cases should guide Stage 5 interaction analysis and Stage 6 feature engineering.
