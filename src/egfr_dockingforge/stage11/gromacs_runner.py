from __future__ import annotations

import json
import os
import subprocess
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import pandas as pd

from egfr_dockingforge.common.io import write_table
from egfr_dockingforge.stage11.forcefield_config import gromacs_version


# Shared equilibration phases (run once per system).
EQUIL_PHASES = [
    ("minimization", "EM"),
    ("nvt_equilibration", "NVT"),
    ("npt_equilibration", "NPT"),
]
PRODUCTION_PHASE = ("production_quick", "NPT")

# Back-compat: some callers/tests import PHASES.
PHASES = [*EQUIL_PHASES, PRODUCTION_PHASE]


def _num_replicates(sys: dict[str, Any], config: dict[str, Any]) -> int:
    """Replicate count for a system: finalists get finalist_replicates, others
    get quick_replicates. Independent replicates use distinct velocity seeds."""
    md = config["md"]
    if bool(sys.get("selected_for_replicate_md", False)):
        return max(1, int(md.get("finalist_replicates", 3)))
    return max(1, int(md.get("quick_replicates", 1)))


def _write_replicate_production_mdp(base_mdp: Path, work: Path, replicate_index: int) -> Path:
    """Create a production .mdp for an independent replicate: regenerate
    velocities from a distinct seed instead of continuing the NPT velocities,
    so replicates are statistically independent."""
    if replicate_index == 1:
        return base_mdp
    text = base_mdp.read_text(encoding="utf-8")
    seed = 20260708 + 1000 * replicate_index
    lines = []
    for line in text.splitlines():
        key = line.split("=")[0].strip().lower()
        if key in {"gen_vel", "gen-vel", "continuation", "gen_seed", "gen-seed"}:
            continue
        lines.append(line)
    lines += [f"gen_vel = yes", f"gen_temp = {int(float(_ref_t(text)))}", f"gen_seed = {seed}", "continuation = no"]
    out = work / f"production_quick_rep{replicate_index:02d}.mdp"
    out.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return out


def _ref_t(mdp_text: str) -> str:
    for line in mdp_text.splitlines():
        if line.split("=")[0].strip().lower() in {"ref_t", "ref-t"}:
            return line.split("=", 1)[1].strip().split()[0]
    return "300"


def _run(command: list[str], log_path: Path) -> tuple[int, str, str]:
    start = datetime.now(timezone.utc)
    result = subprocess.run(command, text=True, capture_output=True)
    end = datetime.now(timezone.utc)
    text = f"$ {' '.join(command)}\n\nSTDOUT:\n{result.stdout}\n\nSTDERR:\n{result.stderr}\n"
    log_path.parent.mkdir(parents=True, exist_ok=True)
    log_path.write_text(text, encoding="utf-8")
    return result.returncode, start.isoformat(), end.isoformat()


def _phase_input(work: Path, sys: dict[str, Any], phase: str) -> str:
    if phase == "minimization":
        return sys["ionized_structure_file"]
    if phase == "nvt_equilibration":
        return str(work / "minimization.gro")
    if phase == "npt_equilibration":
        return str(work / "nvt_equilibration.gro")
    if phase == "production_quick":
        return str(work / "npt_equilibration.gro")
    raise ValueError(phase)


def _phase_extra_grompp(work: Path, phase: str) -> list[str]:
    if phase == "nvt_equilibration":
        return ["-r", str(work / "minimization.gro")]
    if phase == "npt_equilibration":
        return ["-t", str(work / "nvt_equilibration.cpt"), "-r", str(work / "nvt_equilibration.gro")]
    if phase == "production_quick":
        return ["-t", str(work / "npt_equilibration.cpt")]
    return []


def _phase_complete(work: Path, phase: str) -> bool:
    if phase == "production_quick":
        return (work / f"{phase}.xtc").exists() and (work / f"{phase}.gro").exists()
    return (work / f"{phase}.gro").exists()


def _record_run(sys, phase, ensemble, replicate_id, work, deffnm, mdp, config, version, status, exit_code, start_time, end_time, error, input_structure) -> dict[str, Any]:
    gmx = config["forcefield"]["gromacs_executable"]
    tpr = work / f"{deffnm}.tpr"
    grompp_line = f"{gmx} grompp -f {mdp} -c {input_structure} -p {sys.get('topology_file','')} -o {tpr} -maxwarn 0"
    return {
        "md_run_id": f"run_{sys['md_system_id']}_{phase}_{replicate_id}",
        "md_system_id": sys["md_system_id"],
        "md_candidate_id": sys["md_candidate_id"],
        "replicate_id": replicate_id,
        "md_phase": phase,
        "ensemble": ensemble,
        "input_structure": input_structure,
        "topology_file": sys["topology_file"],
        "mdp_file": str(mdp),
        "tpr_file": str(tpr),
        "output_structure": str(work / f"{deffnm}.gro"),
        "trajectory_file": str(work / f"{deffnm}.xtc"),
        "energy_file": str(work / f"{deffnm}.edr"),
        "log_file": str(work / f"{deffnm}.log"),
        "checkpoint_file": str(work / f"{deffnm}.cpt"),
        "gromacs_version": version,
        "command_line_grompp": grompp_line,
        "command_line_mdrun": f"{gmx} mdrun -deffnm {work / deffnm}",
        "start_time": start_time,
        "end_time": end_time,
        "runtime_seconds": 0.0,
        "gpu_used_flag": config["md"]["gpu"] != "cpu",
        "exit_code": exit_code,
        "run_status": status,
        "error_message": error,
        "warnings_json": json.dumps([] if status == "complete" else [status]),
    }


def _run_phase(gmx, work, deffnm, mdp, topology, input_structure, extra_grompp, mdrun_extra=None) -> tuple[str, int | None, str, str, str]:
    """grompp + mdrun for one phase/replicate with a given deffnm. Returns
    (status, exit_code, start, end, error).

    Power-outage resilient:
      * If the phase already finished (.gro, and .xtc for production), skip it.
      * If a checkpoint (.cpt) exists but the run did not finish, resume the
        mdrun from the checkpoint with `-cpi` instead of restarting.
      * grompp is skipped when the .tpr already exists (it is deterministic).

    ``mdrun_extra`` (list) appends thread/pinning flags so multiple runs can
    share the GPU concurrently without oversubscribing cores.
    """
    tpr = work / f"{deffnm}.tpr"
    gro = work / f"{deffnm}.gro"
    xtc = work / f"{deffnm}.xtc"
    cpt = work / f"{deffnm}.cpt"
    is_prod = deffnm.startswith("production")
    if gro.exists() and (not is_prod or xtc.exists()):
        return "complete", 0, "", "", ""
    start_time = ""
    end_time = ""
    if not tpr.exists():
        grompp = [gmx, "grompp", "-f", str(mdp), "-c", input_structure, "-p", topology, "-o", str(tpr), "-maxwarn", "0"] + extra_grompp
        grompp_code, start_time, end_time = _run(grompp, work / f"{deffnm}.grompp.exec.log")
        if grompp_code != 0:
            return "failed_grompp", grompp_code, start_time, end_time, f"grompp failed; see {work / f'{deffnm}.grompp.exec.log'}"
    mdrun = [gmx, "mdrun", "-deffnm", str(work / deffnm)]
    if mdrun_extra:
        mdrun += list(mdrun_extra)
    if cpt.exists():
        # Resume an interrupted run from its last checkpoint.
        mdrun += ["-cpi", str(cpt)]
    mdrun_code, md_start, md_end = _run(mdrun, work / f"{deffnm}.mdrun.exec.log")
    if mdrun_code != 0:
        return "failed_mdrun", mdrun_code, start_time or md_start, md_end, f"mdrun failed; see {work / f'{deffnm}.mdrun.exec.log'}"
    return "complete", 0, start_time or md_start, md_end, ""


def _pin_flags(config: dict[str, Any], slot: int) -> list[str]:
    """Per-run OpenMP thread + core-pin flags so N concurrent mdruns on one GPU
    do not oversubscribe CPU cores. Cores are partitioned into `max_concurrent`
    contiguous blocks; slot i gets its own block. GPU is shared (small systems
    under-utilize the card, so co-scheduling raises aggregate throughput)."""
    md = config["md"]
    total = int(md.get("cpu_cores", os.cpu_count() or 12))
    conc = max(1, int(md.get("max_concurrent_runs", 1)))
    ntomp = max(2, total // conc)
    flags = ["-ntomp", str(ntomp)]
    if bool(md.get("pin_threads", conc > 1)):
        flags += ["-pin", "on", "-pinoffset", str((slot % conc) * ntomp), "-pinstride", "1"]
    # Force this run onto the single GPU explicitly when co-scheduling.
    if md.get("gpu") not in ("cpu",):
        flags += ["-gpu_id", str(md.get("gpu_id", "0"))]
    return flags


def _equilibrate_system(sys, config, paths, version) -> tuple[list[dict], bool]:
    """Run the (internally serial) min->NVT->NPT chain for one system. Returns
    (recorded rows, equilibration_ok)."""
    gmx = config["forcefield"]["gromacs_executable"]
    work = paths["md_root"] / sys["md_system_id"]
    topology = sys["topology_file"]
    blocked = sys["build_status"] != "ready_for_md"
    rows = []
    if blocked:
        phase, ensemble = EQUIL_PHASES[0]
        rows.append(_record_run(sys, phase, ensemble, "rep01", work, phase, work / f"{phase}.mdp", config, version, "blocked_missing_system_build", None, "", "", "system build blocked", ""))
        return rows, False
    for phase, ensemble in EQUIL_PHASES:
        mdp = work / f"{phase}.mdp"
        input_structure = _phase_input(work, sys, phase)
        status, code, st, en, err = _run_phase(gmx, work, phase, mdp, topology, input_structure, _phase_extra_grompp(work, phase))
        rows.append(_record_run(sys, phase, ensemble, "rep01", work, phase, mdp, config, version, status, code, st, en, err, input_structure))
        if status != "complete":
            return rows, False
    return rows, True


def make_gromacs_runs(systems: pd.DataFrame, config: dict[str, Any], paths: dict[str, Path]) -> pd.DataFrame:
    from concurrent.futures import ThreadPoolExecutor, as_completed

    gmx = config["forcefield"]["gromacs_executable"]
    version = gromacs_version(gmx)
    conc = max(1, int(config["md"].get("max_concurrent_runs", 1)))
    # Equilibration (esp. steepest-descent minimization) has tiny per-step GPU
    # work, so high concurrency thrashes the GPU launch queue and is far slower
    # than giving each run more CPU threads. Use a lower equilibration
    # concurrency (default 2) so min/NVT/NPT finish quickly; reserve full
    # concurrency for the long production phase where it gives the ~3x win.
    equil_conc = max(1, int(config["md"].get("max_concurrent_equil", min(2, conc))))
    rows: list[dict] = []

    sys_records = systems.to_dict("records")

    # --- Phase 1: equilibrate systems (each chain internal-serial; systems in
    # parallel up to `equil_conc`). Only systems whose equilibration completes proceed. ---
    equil_ok: dict[str, bool] = {}
    with ThreadPoolExecutor(max_workers=equil_conc) as pool:
        futs = {pool.submit(_equilibrate_system, sys, config, paths, version): sys for sys in sys_records}
        for fut in as_completed(futs):
            sys = futs[fut]
            erows, ok = fut.result()
            rows.extend(erows)
            equil_ok[sys["md_system_id"]] = ok

    # --- Phase 2: all production replicates across all ready systems, run through
    # one bounded pool so up to `conc` share the GPU concurrently. ---
    phase, ensemble = PRODUCTION_PHASE
    jobs = []  # (sys, replicate_id, deffnm, mdp, extra)
    for sys in sys_records:
        if not equil_ok.get(sys["md_system_id"]):
            continue
        work = paths["md_root"] / sys["md_system_id"]
        base_mdp = work / f"{phase}.mdp"
        n_rep = _num_replicates(sys, config)
        for r in range(1, n_rep + 1):
            deffnm = phase if r == 1 else f"{phase}_rep{r:02d}"
            mdp = _write_replicate_production_mdp(base_mdp, work, r)
            extra = ["-t", str(work / "npt_equilibration.cpt")] if r == 1 else ["-r", str(work / "npt_equilibration.gro")]
            jobs.append((sys, f"rep{r:02d}", deffnm, mdp, extra))

    def _run_production(job, slot):
        sys, replicate_id, deffnm, mdp, extra = job
        work = paths["md_root"] / sys["md_system_id"]
        input_structure = str(work / "npt_equilibration.gro")
        status, code, st, en, err = _run_phase(
            gmx, work, deffnm, mdp, sys["topology_file"], input_structure, extra,
            mdrun_extra=_pin_flags(config, slot),
        )
        return _record_run(sys, phase, ensemble, replicate_id, work, deffnm, mdp, config, version, status, code, st, en, err, input_structure)

    with ThreadPoolExecutor(max_workers=conc) as pool:
        futs = {pool.submit(_run_production, job, i): job for i, job in enumerate(jobs)}
        for fut in as_completed(futs):
            rows.append(fut.result())

    out = pd.DataFrame(rows)
    write_table(paths["processed"] / "md_runs.parquet", out)
    write_table(paths["processed"] / "md_runs.csv", out)
    return out
