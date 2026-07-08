# Stage 3 Pose Labeling Rules

Stage 3 labels are provisional and RMSD/sanity based:

- `strict_native_like`: RMSD <= 2.0 A and sanity passes.
- `relaxed_native_like`: RMSD <= 3.0 A and sanity passes.
- `sampled_not_ranked`: native-like pose exists in top-N but top1 is not native-like.
- `ranking_failure`: sampling succeeds but top-ranked pose is wrong.
- `sampling_failure`: no relaxed native-like pose in top-N.
- `invalid_pose`: severe sanity failure.
- `failed_run`: docking/prep/run failed.
- `pending_review`: insufficient data.

Interaction recovery remains `pending_stage5`.

