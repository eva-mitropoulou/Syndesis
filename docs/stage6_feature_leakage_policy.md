# Stage 6 Feature Leakage Policy

Forbidden deployment features include RMSD fields, native-reference identifiers, strict/relaxed native-like flags, final pose labels, Stage 3 pose labels, native IFP similarity, and native key-interaction recovery metrics.

The leakage audit writes one row per candidate feature and marks the action as `train`, `metadata_only`, or `drop_for_training`. Stage 6 training fails if any forbidden feature reaches the model feature schema.
