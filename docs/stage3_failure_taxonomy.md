# Stage 3 Failure Taxonomy

- `engine_unavailable`: requested docking executable/container is unavailable.
- `preparation_failure`: receptor or ligand docking input could not be prepared.
- `sampling_failure`: no near-native pose appears in top-N.
- `ranking_failure`: near-native pose exists but top pose is wrong.
- `invalid_pose`: physical sanity checks fail.
- `state_mismatch`: expected cross-state incompatibility.
- `pending`: task has not been run yet.

