# Stage 6 Pose Model Protocol

Stage 6 trains pose reranking and calibrated pose-confidence models from Stage 3 RMSD/sanity labels, Stage 4 docking/GNINA scores, and Stage 5 ProLIF interaction features.

The deployed feature matrix excludes native-pose comparison fields. RMSD, native IFP similarity, native key-interaction recovery, and final labels are used only to build labels, diagnostics, and retrospective evaluation.

Primary training uses grouped learning-to-rank by `docking_task_id`. Groups without relevance variation are retained for diagnostics but excluded from ranker fitting.
