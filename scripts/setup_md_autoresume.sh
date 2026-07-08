#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# Install auto-resume for the finalist MD driver so it survives power outages
# and abrupt reboots on this WSL2 host.
#
# WSL2 does not run systemd by default and does not honour cron @reboot unless a
# service starts cron. This script therefore installs TWO complementary hooks:
#
#   1. A cron @reboot entry (used if cron is running).
#   2. A ~/.bashrc guard that (re)launches the driver on the first interactive
#      login after a reboot -- practical on WSL where you open a shell to use it.
#
# Both call the SAME idempotent driver (run_finalist_md.sh), which no-ops if the
# work is already complete or another instance holds the flock. Safe to run more
# than once.
#
# Usage:  bash scripts/setup_md_autoresume.sh install
#         bash scripts/setup_md_autoresume.sh status
#         bash scripts/setup_md_autoresume.sh start     # launch now in background
# ---------------------------------------------------------------------------
set -u
PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
DRIVER="${PROJECT_DIR}/scripts/run_finalist_md.sh"
LOGDIR="${PROJECT_DIR}/data/processed/stage11/md_work/logs"
GUARD_MARKER="# >>> egfr finalist MD autoresume >>>"

install_cron() {
  # Add an @reboot entry if not already present.
  local entry="@reboot /usr/bin/env bash ${DRIVER} >> ${LOGDIR}/autoresume_cron.log 2>&1"
  ( crontab -l 2>/dev/null | grep -vF "$DRIVER"; echo "$entry" ) | crontab - 2>/dev/null \
    && echo "cron @reboot installed" || echo "cron not available (bashrc guard will cover it)"
}

install_bashrc_guard() {
  local rc="$HOME/.bashrc"
  if grep -qF "$GUARD_MARKER" "$rc" 2>/dev/null; then
    echo ".bashrc guard already present"
    return
  fi
  cat >> "$rc" <<EOF

${GUARD_MARKER}
# On first interactive login, resume the finalist MD if it hasn't finished and
# isn't already running. Non-blocking: launches detached with nohup.
if [ -x "${DRIVER}" ] && [ ! -f "${LOGDIR}/finalist_md_driver.DONE" ]; then
  if ! pgrep -f run_finalist_md.sh >/dev/null 2>&1; then
    nohup bash "${DRIVER}" >> "${LOGDIR}/autoresume_login.log" 2>&1 &
    disown 2>/dev/null || true
  fi
fi
# <<< egfr finalist MD autoresume <<<
EOF
  echo ".bashrc guard installed"
}

case "${1:-status}" in
  install)
    mkdir -p "$LOGDIR"
    chmod +x "$DRIVER" 2>/dev/null || true
    install_cron
    install_bashrc_guard
    echo "Auto-resume installed. The MD will resume on reboot or next login."
    ;;
  start)
    mkdir -p "$LOGDIR"
    chmod +x "$DRIVER"
    nohup bash "$DRIVER" >> "${LOGDIR}/autoresume_manual.log" 2>&1 &
    disown 2>/dev/null || true
    echo "Driver launched in background (pid $!). tail -f ${LOGDIR}/finalist_md_driver.log"
    ;;
  status)
    if pgrep -f run_finalist_md.sh >/dev/null 2>&1; then echo "driver: RUNNING"; else echo "driver: not running"; fi
    [ -f "${LOGDIR}/finalist_md_driver.DONE" ] && echo "state: COMPLETE" || echo "state: incomplete"
    echo "--- last 15 log lines ---"; tail -n 15 "${LOGDIR}/finalist_md_driver.log" 2>/dev/null
    ;;
  *)
    echo "usage: $0 {install|start|status}"; exit 1;;
esac
