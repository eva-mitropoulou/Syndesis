#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# Resilient finalist MD driver for Syndesis (Stage 11).
#
# Runs the full finalist MD pipeline sequentially and is safe to re-run after a
# power outage or reboot: every step is idempotent (GROMACS phases skip when
# already complete; interrupted mdrun resumes from its .cpt checkpoint via -cpi).
#
# Usage:
#   scripts/run_finalist_md.sh                 # run / resume
#   tail -f data/processed/stage11/md_work/logs/finalist_md_driver.log
#
# Designed to be launched by cron @reboot and/or systemd (see the generated
# unit / crontab suggestions in scripts/setup_md_autoresume.sh).
# ---------------------------------------------------------------------------
set -u

# Resolve the project root from this script's location (scripts/ is at repo root),
# so the driver works regardless of where the repo is checked out.
PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
PY="${PROJECT_DIR}/.tools/stage11_amber/bin/python"
CONFIG="configs/stage11_md_stress_test.yaml"
LOGDIR="${PROJECT_DIR}/data/processed/stage11/md_work/logs"
LOG="${LOGDIR}/finalist_md_driver.log"
LOCK="${LOGDIR}/finalist_md_driver.lock"
STATE="${LOGDIR}/finalist_md_driver.state"

cd "$PROJECT_DIR" || exit 1
mkdir -p "$LOGDIR"

# --- single-instance lock (survives reboot; flock auto-releases on process death) ---
exec 9>"$LOCK"
if ! flock -n 9; then
  echo "$(date -Is) another finalist MD driver instance is already running; exiting." >>"$LOG"
  exit 0
fi

log() { echo "$(date -Is) $*" | tee -a "$LOG"; }

export PYTHONPATH=src
export GMX_MAXBACKUP=-1
export PATH="$PATH:/usr/lib/wsl/lib"

# Minimum free space (GB) required before starting a heavy MD step. Each 20 ns
# production replicate needs ~1 GB; we keep a safety margin so an outage never
# happens because the disk filled. Override with EGFR_MD_MIN_FREE_GB.
MIN_FREE_GB="${EGFR_MD_MIN_FREE_GB:-4}"

free_gb() { df -BG --output=avail / | tail -1 | tr -dc '0-9'; }

require_disk() {
  local avail; avail="$(free_gb)"
  if [ "${avail:-0}" -lt "$MIN_FREE_GB" ]; then
    log "ABORT low disk: ${avail}G free < ${MIN_FREE_GB}G required. Free space or expand the WSL disk (see docs/reproducibility_notes.md), then re-run; the driver resumes where it left off."
    exit 3
  fi
}

log "==== finalist MD driver start (pid $$) ===="
log "free disk: $(free_gb)G (min required ${MIN_FREE_GB}G)"

# Each step writes a marker in STATE when done, so re-runs skip finished steps
# quickly. The heavy GROMACS steps are ALSO internally idempotent, so even a
# missing marker only re-checks (does not recompute) completed work.
step_done() { grep -qxF "$1" "$STATE" 2>/dev/null; }
mark_done() { echo "$1" >>"$STATE"; }

run_step() {
  local name="$1"; shift
  if step_done "$name"; then
    log "SKIP  $name (already marked done)"
    return 0
  fi
  log "START $name : $*"
  if "$@" >>"$LOG" 2>&1; then
    mark_done "$name"
    log "OK    $name"
    return 0
  else
    log "FAIL  $name (exit $?). Will retry on next driver run."
    return 1
  fi
}

# 1. Select finalists (controls + top-K accepted analogs). Cheap, always safe.
run_step select_candidates      "$PY" -m syndesis.cli select-md-candidates --config "$CONFIG" || exit 1
# 2. Parameterize finalist ligands (ACPYPE/GAFF2/AM1-BCC; convergence-gated).
run_step parameterize_ligands   "$PY" -m syndesis.cli parameterize-md-ligands --config "$CONFIG" || exit 1
# 3. Build solvated/ionized systems for finalists.
run_step build_systems          "$PY" -m syndesis.cli build-md-systems --config "$CONFIG" || exit 1
require_disk   # guard before the heavy, disk-hungry production phase
# 4. Run minimization + equilibration + 3 production replicates per finalist.
#    Internally resumable (skips complete phases; -cpi resumes interrupted runs).
#    NOT marked done until it fully completes, so an outage mid-production simply
#    re-enters and continues.
run_step run_md                 "$PY" -m syndesis.cli run-md-production --config "$CONFIG" || exit 1
# 5. Analysis (PBC + backbone-superposition RMSD, min-image interactions).
run_step analyze_trajectories   "$PY" -m syndesis.cli analyze-md-trajectories --config "$CONFIG" || exit 1
run_step interaction_persistence "$PY" -m syndesis.cli compute-md-interaction-persistence --config "$CONFIG" || exit 1
run_step score_stability        "$PY" -m syndesis.cli score-md-stability --config "$CONFIG" || exit 1
run_step report_stage11         "$PY" -m syndesis.cli report-stage11 --config "$CONFIG" || exit 1
# 6. Refresh final dossiers with the new MD evidence.
run_step run_stage12            "$PY" -m syndesis.cli run-stage12 --config configs/stage12_candidate_dossiers.yaml || exit 1

# 7. Auto-finalize: regenerate EVERY stage report (5,6,9,10,11,12) from the
#    corrected pipeline so no on-disk artifact contradicts the manuscript. This
#    fires automatically the instant the MD completes (incl. overnight / after an
#    outage-driven resume), so the project reaches a fully-consistent, portfolio-
#    ready state with no manual step. Idempotent; report-only.
if [ -x "${PROJECT_DIR}/scripts/finalize_all_reports.sh" ]; then
  if ! grep -qxF "finalize_reports" "$STATE" 2>/dev/null; then
    log "START finalize_all_reports"
    if bash "${PROJECT_DIR}/scripts/finalize_all_reports.sh" >>"$LOG" 2>&1; then
      mark_done "finalize_reports"; log "OK    finalize_all_reports"
    else
      log "FAIL  finalize_all_reports (exit $?); reports may be partially updated"
    fi
  else
    log "SKIP  finalize_all_reports (already done)"
  fi
fi

log "==== finalist MD driver COMPLETE ===="
touch "${LOGDIR}/finalist_md_driver.DONE"
