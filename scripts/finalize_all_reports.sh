#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# Regenerate EVERY stage report from the corrected pipeline so no on-disk
# artifact contradicts the manuscript. Run after the finalist MD completes.
#
# Idempotent and safe to re-run. Regenerates stages 5,6,9,10 (from existing
# data) and 11,12 (from the finished MD), then rebuilds the top-level summary.
# Re-reports deterministic Stage 9 analog outputs without regenerating molecules.
# ---------------------------------------------------------------------------
set -u
PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
PY="${PROJECT_DIR}/.tools/stage11_amber/bin/python"
LOG="${PROJECT_DIR}/data/processed/stage11/md_work/logs/finalize_reports.log"
cd "$PROJECT_DIR" || exit 1
export PYTHONPATH=src
export PATH="$PATH:/usr/lib/wsl/lib"
log(){ echo "$(date -Is) $*" | tee -a "$LOG"; }

run(){ log "RUN $*"; if "$PY" -m egfr_dockingforge.cli "$@" >>"$LOG" 2>&1; then log "OK  $1"; else log "FAIL $1 (exit $?)"; fi; }

log "==== finalize_all_reports start ===="

# Stage 5 interaction recovery (conserved-core consensus) + report
run compute-interaction-recovery --config configs/stage5_interaction_atlas.yaml
run report-stage5 --config configs/stage5_interaction_atlas.yaml

# Stage 6 leakage-safe pose model report
run build-pose-features   --config configs/stage6_pose_model.yaml
run audit-pose-features   --config configs/stage6_pose_model.yaml
run build-pose-labels     --config configs/stage6_pose_model.yaml
run evaluate-pose-models  --config configs/stage6_pose_model.yaml
run report-stage6         --config configs/stage6_pose_model.yaml

# Stage 9 re-report from existing deterministic analog outputs
run score-analog-acceptance      --config configs/stage9_deterministic_analogs.yaml
run benchmark-analog-strategies  --config configs/stage9_deterministic_analogs.yaml
run report-stage9                --config configs/stage9_deterministic_analogs.yaml

# Stage 10 ablation (consistent acceptance tolerances) + report
run compute-analog-benchmark-metrics --config configs/stage10_ablation_benchmark.yaml
run compute-score-hacking-metrics    --config configs/stage10_ablation_benchmark.yaml
run run-ablation-statistics          --config configs/stage10_ablation_benchmark.yaml
run report-stage10                   --config configs/stage10_ablation_benchmark.yaml

# Stage 11 MD (analysis already done by the MD driver; re-report to be safe)
run analyze-md-trajectories            --config configs/stage11_md_stress_test.yaml
run compute-md-interaction-persistence --config configs/stage11_md_stress_test.yaml
run score-md-stability                 --config configs/stage11_md_stress_test.yaml
run report-stage11                     --config configs/stage11_md_stress_test.yaml

# Stage 12 final dossiers + provenance bundle (consumes corrected MD)
run run-stage12 --config configs/stage12_candidate_dossiers.yaml

log "==== finalize_all_reports COMPLETE ===="
touch "${PROJECT_DIR}/data/processed/stage11/md_work/logs/finalize_reports.DONE"
