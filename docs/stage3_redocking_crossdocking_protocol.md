# Stage 3 Redocking And Cross-docking Protocol

Stage 3 builds an engine-agnostic redocking/cross-docking benchmark. It separates
sampling failure, ranking failure, receptor-state mismatch, physically invalid poses, and
tool/preparation failures.

The default baseline engine is AutoDock Vina. If Vina is not installed, the pipeline
records explicit `engine_unavailable` run metadata and does not fabricate poses.

Stage 3 does not implement GNINA rescoring, ProLIF interaction recovery, pose-confidence
modeling, MD, or agentic analog generation.

