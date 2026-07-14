from __future__ import annotations

import json
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd

from egfr_dockingforge.common.io import write_table


# ---------------------------------------------------------------------------
# Trajectory preprocessing: periodic-boundary correction + protein superposition
# ---------------------------------------------------------------------------
#
# Ligand-pose stability MUST be measured in the protein reference frame after
# removing (a) periodic-boundary jumps and (b) global rigid-body rotation and
# translation of the whole solvated complex. Computing an RMSD on raw lab-frame
# coordinates measures the complex diffusing/tumbling through the water box, not
# ligand motion, and yields spuriously huge values (16-21 A for stable poses).
#
# We use ``gmx trjconv`` for this because MDAnalysis cannot parse the GROMACS
# 2026 TPR (tpx version 138). The corrected trajectory is cached next to the
# production trajectory and reused by the interaction-persistence analysis.


def _run_gmx(command: list[str], stdin_text: str, log_path: Path) -> subprocess.CompletedProcess:
    log_path.parent.mkdir(parents=True, exist_ok=True)
    result = subprocess.run(
        command,
        input=stdin_text,
        text=True,
        capture_output=True,
        env={"GMX_MAXBACKUP": "-1", "PATH": _os_path()},
    )
    log_path.write_text(
        f"$ {' '.join(command)}\n\nSTDIN:\n{stdin_text}\n\nSTDOUT:\n{result.stdout}\n\nSTDERR:\n{result.stderr}\n",
        encoding="utf-8",
    )
    return result


def _os_path() -> str:
    import os

    return os.environ.get("PATH", "/usr/local/gromacs/bin:/usr/bin:/bin")


def pbc_corrected_trajectory(gmx: str, tpr: Path, xtc: Path, work: Path) -> tuple[Path, list[str]]:
    """Return a PBC-molecule-corrected, backbone-fitted, water-stripped trajectory.

    Step 1: ``-pbc mol -center`` centres the protein and makes molecules whole.
    Step 2: ``-fit rot+trans`` fits every frame onto the protein Backbone (so
            ligand motion is measured in the protein frame) and writes ONLY the
            non-Water group (protein + ligand + ions). Dropping bulk water shrinks
            the analysis trajectory ~20x with no effect on ligand-pose metrics,
            which keeps the multi-replicate MD within disk budget.

    A matching water-stripped .gro topology (``<deffnm>_noPBCwater.gro``) is
    written alongside so MDAnalysis atom counts line up with the trajectory.
    Cached: re-analysis reuses the fitted trajectory if newer than the source.
    """
    warnings: list[str] = []
    fitted = work / (xtc.stem + "_pbcfit.xtc")
    stripped_gro = work / (xtc.stem + "_noPBCwater.gro")
    fitted_ok = fitted.exists() and fitted.stat().st_size > 0 and fitted.stat().st_mtime >= xtc.stat().st_mtime
    gro_ok = stripped_gro.exists() and stripped_gro.stat().st_mtime >= xtc.stat().st_mtime
    if fitted_ok and gro_ok:
        return fitted, warnings
    if fitted_ok and not gro_ok:
        # Fitted trajectory cached but its matching stripped topology is missing;
        # regenerate just the topology (single frame) from the full-atom source.
        rc = _run_gmx(
            [gmx, "trjconv", "-s", str(tpr), "-f", str(xtc), "-o", str(stripped_gro), "-dump", "0"],
            "non-Water\n",
            work / "trjconv_stripgro.log",
        )
        if rc.returncode != 0 or not stripped_gro.exists():
            warnings.append("trjconv_stripgro_failed")
            raise RuntimeError(f"gmx trjconv non-Water .gro export failed; see {work / 'trjconv_stripgro.log'}")
        return fitted, warnings
    pbc = work / (xtc.stem + "_pbcmol.xtc")
    # Step 1: make whole + centre on protein, keep all atoms (System).
    step1 = _run_gmx(
        [gmx, "trjconv", "-s", str(tpr), "-f", str(xtc), "-o", str(pbc), "-pbc", "mol", "-center", "-ur", "compact"],
        "Protein\nSystem\n",
        work / "trjconv_pbc.log",
    )
    if step1.returncode != 0 or not pbc.exists():
        warnings.append("trjconv_pbc_failed")
        raise RuntimeError(f"gmx trjconv -pbc mol failed; see {work / 'trjconv_pbc.log'}")
    # Step 2: fit on Backbone, output only non-Water atoms.
    step2 = _run_gmx(
        [gmx, "trjconv", "-s", str(tpr), "-f", str(pbc), "-o", str(fitted), "-fit", "rot+trans"],
        "Backbone\nnon-Water\n",
        work / "trjconv_fit.log",
    )
    pbc.unlink(missing_ok=True)
    if step2.returncode != 0 or not fitted.exists():
        warnings.append("trjconv_fit_failed")
        raise RuntimeError(f"gmx trjconv -fit rot+trans failed; see {work / 'trjconv_fit.log'}")
    # Matching water-stripped topology (single frame) for MDAnalysis. Derive it
    # from the FULL-ATOM source (tpr + original xtc) selecting non-Water, so the
    # index is valid; selecting against the already-stripped fitted xtc would
    # mismatch. Its atom set then equals the water-stripped fitted trajectory.
    stripped_gro = work / (xtc.stem + "_noPBCwater.gro")
    if not (stripped_gro.exists() and stripped_gro.stat().st_mtime >= xtc.stat().st_mtime):
        rc = _run_gmx(
            [gmx, "trjconv", "-s", str(tpr), "-f", str(xtc), "-o", str(stripped_gro), "-dump", "0"],
            "non-Water\n",
            work / "trjconv_stripgro.log",
        )
        if rc.returncode != 0 or not stripped_gro.exists():
            warnings.append("trjconv_stripgro_failed")
            raise RuntimeError(f"gmx trjconv non-Water .gro export failed; see {work / 'trjconv_stripgro.log'}")
    return fitted, warnings


# ---------------------------------------------------------------------------
# Equilibration QC from real GROMACS energy files
# ---------------------------------------------------------------------------


def _edr_series(gmx: str, edr: Path, term: str, work: Path) -> np.ndarray:
    """Extract one energy term time series from an .edr via ``gmx energy``."""
    if not edr.exists():
        return np.array([])
    xvg = work / f"_qc_{term.lower().replace('-', '_').replace('.', '')}.xvg"
    result = _run_gmx(
        [gmx, "energy", "-f", str(edr), "-o", str(xvg)],
        f"{term}\n\n",
        work / "gmx_energy.log",
    )
    if result.returncode != 0 or not xvg.exists():
        return np.array([])
    values = []
    for line in xvg.read_text(encoding="utf-8", errors="ignore").splitlines():
        line = line.strip()
        if not line or line.startswith(("#", "@")):
            continue
        parts = line.split()
        if len(parts) >= 2:
            try:
                values.append(float(parts[1]))
            except ValueError:
                continue
    xvg.unlink(missing_ok=True)
    return np.asarray(values, dtype=float)


def _em_fmax_from_log(log: Path) -> float | None:
    if not log.exists():
        return None
    fmax = None
    for line in log.read_text(encoding="utf-8", errors="ignore").splitlines():
        stripped = line.strip()
        if stripped.startswith("Maximum force"):
            # e.g. "Maximum force     =  8.1e+02 on atom 123"
            try:
                fmax = float(stripped.split("=", 1)[1].split()[0])
            except (IndexError, ValueError):
                continue
    return fmax


def _equilibration_qc(row: dict, gmx: str, config: dict) -> dict:
    work = Path(row["output_structure"]).parent
    target_t = float(config["md"]["temperature_k"])
    nvt_edr = work / "nvt_equilibration.edr"
    npt_edr = work / "npt_equilibration.edr"
    em_edr = work / "minimization.edr"
    em_log = work / "minimization.log"

    nvt_temp = _edr_series(gmx, nvt_edr, "Temperature", work)
    npt_temp = _edr_series(gmx, npt_edr, "Temperature", work)
    npt_pres = _edr_series(gmx, npt_edr, "Pressure", work)
    npt_dens = _edr_series(gmx, npt_edr, "Density", work)
    em_pot = _edr_series(gmx, em_edr, "Potential", work)
    em_fmax = _em_fmax_from_log(em_log)

    warnings: list[str] = []

    def _mean(a: np.ndarray) -> float | None:
        return float(np.mean(a)) if a.size else None

    def _std(a: np.ndarray) -> float | None:
        return float(np.std(a)) if a.size else None

    npt_temp_mean = _mean(npt_temp)
    npt_dens_mean = _mean(npt_dens)
    npt_dens_std = _std(npt_dens)

    # Physically-motivated equilibration acceptance criteria.
    temp_ok = npt_temp_mean is not None and abs(npt_temp_mean - target_t) <= 5.0
    dens_ok = npt_dens_mean is not None and 950.0 <= npt_dens_mean <= 1050.0
    dens_stable = npt_dens_std is not None and npt_dens_std <= 15.0
    em_ok = (em_pot.size > 0 and np.isfinite(em_pot[-1])) and (em_fmax is None or em_fmax <= 5000.0)
    minimization_status = "complete" if em_pot.size > 0 else "not_run"

    if npt_temp_mean is None:
        warnings.append("npt_energy_file_unreadable")
    if npt_temp_mean is not None and not temp_ok:
        warnings.append(f"npt_temperature_off_target:{npt_temp_mean:.1f}K_vs_{target_t:.0f}K")
    if npt_dens_mean is not None and not dens_ok:
        warnings.append(f"npt_density_out_of_range:{npt_dens_mean:.1f}")
    if npt_dens_std is not None and not dens_stable:
        warnings.append(f"npt_density_unstable_std:{npt_dens_std:.1f}")
    if not em_ok:
        warnings.append("minimization_did_not_converge")

    equilibration_pass = bool(temp_ok and dens_ok and dens_stable and em_ok)

    return {
        "minimization_status": minimization_status,
        "final_em_potential_energy": float(em_pot[-1]) if em_pot.size else None,
        "final_em_fmax": em_fmax,
        "nvt_temperature_mean": _mean(nvt_temp),
        "nvt_temperature_std": _std(nvt_temp),
        "npt_temperature_mean": npt_temp_mean,
        "npt_pressure_mean": _mean(npt_pres),
        "npt_pressure_std": _std(npt_pres),
        "npt_density_mean": npt_dens_mean,
        "npt_density_std": npt_dens_std,
        "box_volume_stable_flag": bool(dens_stable),
        "equilibration_pass_flag": equilibration_pass,
        "failure_reason": "" if equilibration_pass else ";".join(warnings) or "equilibration_criteria_not_met",
        "warnings_json": json.dumps(warnings),
    }


# ---------------------------------------------------------------------------
# Production trajectory analysis (ligand pose stability in protein frame)
# ---------------------------------------------------------------------------


def _analyze_completed_run(row: dict, gmx: str) -> dict:
    import MDAnalysis as mda
    from MDAnalysis.analysis import align, rms
    from MDAnalysis.lib.distances import distance_array

    tpr = Path(row["tpr_file"])
    xtc = Path(row["trajectory_file"])
    gro = Path(row.get("output_structure", ""))
    work = xtc.parent
    if not xtc.exists() or not gro.exists():
        raise FileNotFoundError("production trajectory or GRO topology missing")
    if not tpr.exists():
        raise FileNotFoundError("production TPR missing (needed for PBC correction)")

    warnings: list[str] = []
    fitted, fit_warnings = pbc_corrected_trajectory(gmx, tpr, xtc, work)
    warnings.extend(fit_warnings)

    # Topology from the water-stripped .gro (MDAnalysis cannot parse the GROMACS
    # 2026 TPR); its atom set matches the water-stripped fitted trajectory.
    stripped_gro = work / (xtc.stem + "_noPBCwater.gro")
    topology = stripped_gro if stripped_gro.exists() else gro
    universe = mda.Universe(str(topology), str(fitted))
    protein = universe.select_atoms("protein")
    backbone = universe.select_atoms("protein and backbone")
    ligand = universe.select_atoms("resname UNL")
    if len(ligand) == 0:
        ligand = universe.select_atoms("not protein and not resname SOL NA CL HOH WAT")
    if len(protein) == 0 or len(ligand) == 0:
        raise ValueError("protein or ligand atom selection is empty")

    # Reference = first frame of the corrected trajectory.
    universe.trajectory[0]
    ligand_ref = ligand.positions.copy()
    backbone_ref = backbone.positions.copy()
    lig_heavy_ref = ligand.select_atoms("not name H*") or ligand
    box0 = universe.dimensions
    # Pocket-lining protein backbone atoms: backbone atoms of residues within 6 A
    # of any ligand heavy atom in the reference frame. This LOCAL, rigid selection
    # is the correct alignment frame for pose analysis. The trjconv global-Backbone
    # fit is dominated by mobile N/C-terminal tails and flexible loops of the
    # kinase domain (backbone RMSD 4-7 A even for a stable pose), which leaves the
    # binding site sub-optimally aligned and INFLATES the ligand RMSD. Re-aligning
    # per frame on the pocket backbone gives (a) a meaningful pocket RMSD and (b) a
    # ligand RMSD measured in the frame that actually matters for binding.
    ref_pair = distance_array(protein.positions, lig_heavy_ref.positions, box=box0)
    pocket_res_mask = ref_pair.min(axis=1) <= 6.0
    pocket_residues = protein[pocket_res_mask].residues
    pocket_bb = pocket_residues.atoms.select_atoms("backbone")
    if len(pocket_bb) < 4:  # fall back to full backbone if the local set is too small
        pocket_bb = backbone
        warnings.append("pocket_backbone_too_small_used_full_backbone")
    pocket_bb_ref = pocket_bb.positions.copy()
    pocket_bb_ref_c = pocket_bb_ref - pocket_bb_ref.mean(axis=0)

    # Pocket heavy atoms (for the in-pocket contact test) = pocket residues' heavy atoms.
    pocket_atoms = pocket_residues.atoms.select_atoms("not name H*")
    ligand_com_ref = ligand.center_of_mass()

    ligand_rmsd = []
    pocket_rmsd = []
    backbone_rmsd = []
    com_drift = []
    rg = []
    inside = []
    times = []
    lig_heavy = ligand.select_atoms("not name H*") or ligand
    lig_heavy_ref_pos = lig_heavy.positions.copy()
    for ts in universe.trajectory:
        # Per-frame least-squares superposition on the POCKET backbone (Kabsch),
        # then apply that same rotation/translation to the ligand so its RMSD is
        # measured in the locally-aligned binding-site frame.
        mob = pocket_bb.positions
        mob_com = mob.mean(axis=0)
        R, _ = align.rotation_matrix(mob - mob_com, pocket_bb_ref_c)
        # pocket backbone RMSD after this optimal local alignment
        aligned_bb = (mob - mob_com) @ R.T
        pocket_rmsd.append(float(np.sqrt(np.mean(np.sum((aligned_bb - pocket_bb_ref_c) ** 2, axis=1)))))
        # ligand heavy-atom RMSD in the same local frame
        lig = lig_heavy.positions
        aligned_lig = (lig - mob_com) @ R.T
        aligned_lig_ref = lig_heavy_ref_pos - pocket_bb_ref.mean(axis=0)
        ligand_rmsd.append(float(np.sqrt(np.mean(np.sum((aligned_lig - aligned_lig_ref) ** 2, axis=1)))))
        # global backbone RMSD (trjconv-fitted frame) kept for reference/QC
        backbone_rmsd.append(float(rms.rmsd(backbone.positions, backbone_ref, center=False, superposition=False)))
        com_drift.append(float(np.linalg.norm(ligand.center_of_mass() - ligand_com_ref)))
        rg.append(float(ligand.radius_of_gyration()))
        if len(pocket_atoms):
            d = distance_array(pocket_atoms.positions, lig_heavy.positions, box=ts.dimensions)
            inside.append(bool(d.min() <= 4.5))
        else:
            inside.append(True)
        times.append(float(ts.time))

    analyzed_ns = (max(times) - min(times)) / 1000.0 if times else 0.0
    ligand_rmsd = np.asarray(ligand_rmsd)
    pocket_rmsd = np.asarray(pocket_rmsd)
    return {
        "analyzed_ns": analyzed_ns,
        "ligand_rmsd_median_angstrom": float(np.median(ligand_rmsd)),
        "ligand_rmsd_p95_angstrom": float(np.percentile(ligand_rmsd, 95)),
        "ligand_rmsd_max_angstrom": float(np.max(ligand_rmsd)),
        "pocket_rmsd_median_angstrom": float(np.median(pocket_rmsd)),
        "pocket_rmsd_p95_angstrom": float(np.percentile(pocket_rmsd, 95)),
        "protein_backbone_rmsd_median_angstrom": float(np.median(backbone_rmsd)),
        "ligand_com_drift_median_angstrom": float(np.median(com_drift)),
        "ligand_com_drift_max_angstrom": float(np.max(com_drift)),
        "ligand_rg_median_angstrom": float(np.median(rg)),
        "fraction_frames_inside_pocket": float(np.mean(inside)),
        "ligand_left_pocket_flag": bool(np.mean(inside) < 0.9),
        "ligand_flip_flag": None,
        "trajectory_analysis_status": "complete",
        "warnings_json": json.dumps(warnings),
    }


def analyze_trajectories(runs: pd.DataFrame, paths: dict[str, Path], config: dict) -> tuple[pd.DataFrame, pd.DataFrame]:
    gmx = config["forcefield"]["gromacs_executable"]
    prod = runs[runs["md_phase"].str.startswith("production")]
    rows = []
    qc_rows = []
    for row in prod.to_dict("records"):
        blocked = row["run_status"] != "complete"
        analysis = None
        warnings = ["production_trajectory_missing"]
        status = "failed_missing_trajectory"
        if not blocked:
            try:
                analysis = _analyze_completed_run(row, gmx)
                warnings = json.loads(analysis["warnings_json"])
                status = analysis["trajectory_analysis_status"]
            except Exception as exc:
                warnings = [f"trajectory_analysis_failed:{exc}"]
                status = "failed_analysis_error"
        rows.append(
            {
                "md_candidate_id": row["md_candidate_id"],
                "md_system_id": row["md_system_id"],
                "replicate_id": row["replicate_id"],
                "production_type": row["md_phase"],
                "trajectory_file": row["trajectory_file"],
                "analyzed_ns": 0.0 if analysis is None else analysis["analyzed_ns"],
                "ligand_rmsd_median_angstrom": None if analysis is None else analysis["ligand_rmsd_median_angstrom"],
                "ligand_rmsd_p95_angstrom": None if analysis is None else analysis["ligand_rmsd_p95_angstrom"],
                "ligand_rmsd_max_angstrom": None if analysis is None else analysis["ligand_rmsd_max_angstrom"],
                "pocket_rmsd_median_angstrom": None if analysis is None else analysis["pocket_rmsd_median_angstrom"],
                "pocket_rmsd_p95_angstrom": None if analysis is None else analysis["pocket_rmsd_p95_angstrom"],
                "protein_backbone_rmsd_median_angstrom": None if analysis is None else analysis["protein_backbone_rmsd_median_angstrom"],
                "ligand_com_drift_median_angstrom": None if analysis is None else analysis["ligand_com_drift_median_angstrom"],
                "ligand_com_drift_max_angstrom": None if analysis is None else analysis["ligand_com_drift_max_angstrom"],
                "ligand_rg_median_angstrom": None if analysis is None else analysis["ligand_rg_median_angstrom"],
                "fraction_frames_inside_pocket": 0.0 if analysis is None else analysis["fraction_frames_inside_pocket"],
                "ligand_left_pocket_flag": None if analysis is None else analysis["ligand_left_pocket_flag"],
                "ligand_flip_flag": None if analysis is None else analysis["ligand_flip_flag"],
                "trajectory_analysis_status": status,
                "warnings_json": json.dumps(warnings),
            }
        )
        if blocked:
            qc = {
                "minimization_status": "not_run",
                "final_em_potential_energy": None,
                "final_em_fmax": None,
                "nvt_temperature_mean": None,
                "nvt_temperature_std": None,
                "npt_temperature_mean": None,
                "npt_pressure_mean": None,
                "npt_pressure_std": None,
                "npt_density_mean": None,
                "npt_density_std": None,
                "box_volume_stable_flag": False,
                "equilibration_pass_flag": False,
                "failure_reason": "production_not_started_due_to_setup_block",
                "warnings_json": json.dumps(["parameterization_or_system_build_failed"]),
            }
        else:
            try:
                qc = _equilibration_qc(row, gmx, config)
            except Exception as exc:
                qc = {
                    "minimization_status": "unknown",
                    "final_em_potential_energy": None,
                    "final_em_fmax": None,
                    "nvt_temperature_mean": None,
                    "nvt_temperature_std": None,
                    "npt_temperature_mean": None,
                    "npt_pressure_mean": None,
                    "npt_pressure_std": None,
                    "npt_density_mean": None,
                    "npt_density_std": None,
                    "box_volume_stable_flag": False,
                    "equilibration_pass_flag": False,
                    "failure_reason": f"equilibration_qc_failed:{exc}",
                    "warnings_json": json.dumps([f"equilibration_qc_failed:{exc}"]),
                }
        qc_rows.append(
            {
                "md_system_id": row["md_system_id"],
                "md_candidate_id": row["md_candidate_id"],
                "replicate_id": row["replicate_id"],
                **qc,
            }
        )
    metrics = pd.DataFrame(rows)
    qc = pd.DataFrame(qc_rows)
    write_table(paths["processed"] / "md_metrics.parquet", metrics)
    write_table(paths["processed"] / "md_metrics.csv", metrics)
    write_table(paths["processed"] / "equilibration_qc.parquet", qc)
    write_table(paths["processed"] / "equilibration_qc.csv", qc)
    return metrics, qc
