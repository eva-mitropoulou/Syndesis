from __future__ import annotations

import json
import os
import shutil
import subprocess
from pathlib import Path
from typing import Any

import pandas as pd
from rdkit import Chem

from egfr_dockingforge.common.io import resolve_path, write_table
from egfr_dockingforge.stage11.forcefield_config import acpype_version_from_prefix, ambertools_version_from_prefix, gromacs_version


def _tool_available(value: str) -> bool:
    path = resolve_path(value)
    return path.exists() or shutil.which(value) is not None


def _tool_path(value: str) -> str:
    path = resolve_path(value)
    return str(path) if path.exists() else value


def _stage7_ligand_path(molecule_id: str) -> Path:
    return resolve_path(f"data/processed/stage7/prepared_ligands/prep_{molecule_id}.sdf")


def _formal_charge(sdf_path: Path) -> int:
    mol = Chem.SDMolSupplier(str(sdf_path), removeHs=False)[0]
    if mol is None:
        raise ValueError(f"Could not parse ligand SDF: {sdf_path}")
    return int(sum(atom.GetFormalCharge() for atom in mol.GetAtoms()))


def _run(cmd: list[str], cwd: Path, env: dict[str, str], timeout: int = 10800) -> tuple[int, str]:
    completed = subprocess.run(cmd, cwd=cwd, env=env, capture_output=True, text=True, timeout=timeout, check=False)
    return completed.returncode, (completed.stdout + "\n" + completed.stderr).strip()


def _collect_acpype_warnings(out_dir: Path) -> list[str]:
    warnings: list[str] = []
    log = out_dir / "acpype.log"
    if log.exists():
        text = log.read_text(encoding="utf-8", errors="replace")
        if "ERROR: Antechamber failed" in text:
            warnings.append("acpype_log_antechamber_failed")
        if "ERROR: Tleap failed" in text:
            warnings.append("acpype_log_tleap_failed")
        if "Fatal Error" in text and "sqm" in text:
            warnings.append("acpype_log_sqm_fatal_error")
        if "charge to be balanced" in text:
            warnings.append("acpype_log_charge_balanced_after_parameterization")
        if "ERROR:" in text and not any(w.startswith("acpype_log_") for w in warnings):
            warnings.append("acpype_log_contains_errors_review_log")
    sqm = out_dir / "sqm.out"
    if sqm.exists():
        text = sqm.read_text(encoding="utf-8", errors="replace")
        if "Unable to achieve self consistency" in text:
            warnings.append("sqm_unable_to_achieve_self_consistency")
        if "No convergence in SCF" in text:
            warnings.append("sqm_no_scf_convergence")
        if "QMMM: ERROR!" in text:
            warnings.append("sqm_qmmm_error")
    return warnings


def _amber_env(config: dict[str, Any]) -> dict[str, str]:
    env = os.environ.copy()
    prefix = resolve_path(config["forcefield"].get("ambertools_prefix", ""))
    if prefix.exists():
        env["PATH"] = str(prefix / "bin") + os.pathsep + env.get("PATH", "")
        env["AMBERHOME"] = str(prefix)
    return env


def _parameterize_amber_gaff2(row: dict, config: dict[str, Any], paths: dict[str, Path]) -> dict:
    ff = config["forcefield"]
    work = paths["md_root"] / "ligand_parameters" / row["molecule_id"]
    work.mkdir(parents=True, exist_ok=True)
    input_sdf = resolve_path(row["ligand_file"])
    warnings: list[str] = []
    if not input_sdf.exists():
        repaired = _stage7_ligand_path(row["molecule_id"])
        if repaired.exists():
            warnings.append(f"stage11_manifest_ligand_path_missing_used_stage7_prepared:{row['ligand_file']}")
            input_sdf = repaired
    if not input_sdf.exists():
        raise FileNotFoundError(f"Missing prepared ligand file for {row['molecule_id']}: {row['ligand_file']}")
    charge = _formal_charge(input_sdf)
    mol2 = work / f"{row['molecule_id']}.mol2"
    env = _amber_env(config)
    obabel = _tool_path(ff["openbabel_executable"])
    basename = row["molecule_id"]
    out_dir = work / f"{basename}.acpype"
    itp = out_dir / f"{basename}_GMX.itp"
    gro = out_dir / f"{basename}_GMX.gro"
    frcmod = out_dir / f"{basename}_AC.frcmod"
    log = out_dir / "acpype.log"

    # Idempotent / resumable: AM1-BCC charge derivation is the slow step
    # (~40 min/ligand). If valid ACPYPE outputs already exist and are newer than
    # the input SDF, reuse them instead of recomputing (also survives outages).
    cached = (
        itp.exists() and gro.exists()
        and itp.stat().st_size > 0 and gro.stat().st_size > 0
        and itp.stat().st_mtime >= input_sdf.stat().st_mtime
    )
    if not cached:
        code, text = _run([obabel, "-isdf", str(input_sdf), "-omol2", "-O", str(mol2)], work, env)
        if code != 0 or not mol2.exists():
            raise RuntimeError(f"OpenBabel MOL2 conversion failed for {row['molecule_id']}: {text[-1000:]}")
        acpype = _tool_path(ff["acpype_executable"])
        code, text = _run(
            [acpype, "-i", str(mol2), "-b", basename, "-c", "bcc", "-n", str(charge), "-a", "gaff2", "-o", "gmx", "-f"],
            work,
            env,
        )
        if code != 0 or not itp.exists() or not gro.exists():
            raise RuntimeError(f"ACPYPE GAFF2 parameterization failed for {row['molecule_id']}: {text[-1500:]}")
    else:
        text = "acpype_outputs_reused_from_cache"
        warnings.append("ligand_parameters_reused_from_existing_acpype_outputs")
    if "WARNING" in text:
        warnings.append("acpype_reported_warning_review_log")
    warnings.extend(_collect_acpype_warnings(out_dir))
    return {
        "input_sdf": str(input_sdf),
        "net_charge": charge,
        "itp": str(itp),
        "gro": str(gro),
        "frcmod": str(frcmod),
        "top": str(out_dir / f"{basename}_GMX.top"),
        "log": str(log),
        "warnings": warnings,
    }


def parameterize_ligands(candidates: pd.DataFrame, config: dict[str, Any], paths: dict[str, Path]) -> pd.DataFrame:
    # By default only parameterize MD finalists: AM1-BCC charge derivation is the
    # slow serial step (~40 min/ligand) and non-finalist analogs are never built
    # into MD systems, so parameterizing them wastes hours and idles the GPU.
    # Set selection.parameterize_all_candidates: true to override.
    if (
        "is_md_finalist" in candidates.columns
        and not bool(config.get("selection", {}).get("parameterize_all_candidates", False))
    ):
        finalists = candidates[candidates["is_md_finalist"].fillna(False).astype(bool)]
        if not finalists.empty:
            candidates = finalists
    backend = config["forcefield"].get("ligand_parameterization_backend", "amber_gaff2_acpype")
    if backend == "amber_gaff2_acpype":
        return parameterize_ligands_amber_gaff2(candidates, config, paths)
    return parameterize_ligands_cgenff(candidates, config, paths)


def parameterize_ligands_amber_gaff2(candidates: pd.DataFrame, config: dict[str, Any], paths: dict[str, Path]) -> pd.DataFrame:
    ff = config["forcefield"]
    warnings_stack = []
    if not ff["protein_ff"].lower().startswith("amber"):
        raise ValueError("GAFF2 ligand parameters cannot be combined with non-AMBER protein force field unless an explicit override is implemented.")
    rows = []
    report_rows = []
    for row in candidates.to_dict("records"):
        warnings: list[str] = []
        status = "failed_parameterization"
        reason = ""
        result: dict[str, Any] = {}
        try:
            result = _parameterize_amber_gaff2(row, config, paths)
            warnings.extend(result["warnings"])
            # Gate on AM1-BCC/SQM charge convergence: a ligand whose semi-empirical
            # charge calculation did not converge has unreliable partial charges
            # and must NOT silently proceed to production MD.
            fatal_markers = {
                "acpype_log_sqm_fatal_error",
                "sqm_unable_to_achieve_self_consistency",
                "sqm_no_scf_convergence",
                "sqm_qmmm_error",
                "acpype_log_antechamber_failed",
                "acpype_log_tleap_failed",
            }
            hit = sorted(fatal_markers.intersection(warnings))
            if hit:
                status = "failed_charge_convergence"
                reason = f"AM1-BCC/SQM parameterization not reliable: {', '.join(hit)}"
            else:
                status = "ready"
        except Exception as exc:
            reason = str(exc)
            warnings.append("amber_gaff2_acpype_parameterization_failed")
        rows.append(
            {
                "md_candidate_id": row["md_candidate_id"],
                "molecule_id": row["molecule_id"],
                "ligand_smiles": row["standard_smiles"],
                "ligand_input_file": result.get("input_sdf", row["ligand_file"]),
                "ligand_protonation_state_id": "stage7_prepared_state",
                "ligand_tautomer_state_id": "stage7_prepared_state",
                "cgenff_str_file": "",
                "cgenff_version": "",
                "cgenff_penalty_max": None,
                "cgenff_penalty_mean": None,
                "cgenff_penalty_status": "not_applicable_amber_gaff2_backend",
                "charmm2gmx_status": "not_applicable_amber_gaff2_backend",
                "ligand_itp_file": result.get("itp", ""),
                "ligand_prm_file": result.get("frcmod", ""),
                "ligand_gro_file": result.get("gro", ""),
                "parameterization_status": status,
                "rejection_reason": reason,
                "warnings_json": json.dumps(warnings),
            }
        )
        report_rows.append(
            {
                "ligand_id": row["md_candidate_id"],
                "molecule_id": row["molecule_id"],
                "backend": "amber_gaff2_acpype",
                "input_ligand_file": result.get("input_sdf", row["ligand_file"]),
                "output_itp": result.get("itp", ""),
                "output_gro": result.get("gro", ""),
                "net_charge": result.get("net_charge", None),
                "charge_model": ff.get("ligand_charge_model", "am1bcc"),
                "ligand_forcefield": ff["ligand_ff"],
                "protein_forcefield": ff["protein_ff"],
                "water_model": ff["water_model"],
                "ambertools_version": ambertools_version_from_prefix(ff.get("ambertools_prefix", "")),
                "acpype_version": acpype_version_from_prefix(ff.get("ambertools_prefix", "")),
                "gromacs_version": gromacs_version(ff["gromacs_executable"]),
                "parameterization_status": status,
                "warnings_json": json.dumps(warnings),
            }
        )
    out = pd.DataFrame(rows)
    report = pd.DataFrame(report_rows)
    write_table(paths["processed"] / "ligand_parameterization.parquet", out)
    write_table(paths["processed"] / "ligand_parameterization.csv", out)
    write_table(paths["processed"] / "ligand_parameterization_report.parquet", report)
    write_table(paths["processed"] / "ligand_parameterization_report.csv", report)
    return out


def parameterize_ligands_cgenff(candidates: pd.DataFrame, config: dict[str, Any], paths: dict[str, Path]) -> pd.DataFrame:
    cgenff_found = _tool_available(config["forcefield"]["cgenff_executable"])
    converter_found = _tool_available(config["forcefield"]["charmm2gmx_executable"])
    rows = []
    for row in candidates.to_dict("records"):
        str_file = paths["user_cgenff_str_dir"] / f"{row['molecule_id']}.str"
        has_str = str_file.exists()
        status = "ready" if has_str and converter_found else "failed_parameterization"
        reason = ""
        warnings = []
        if not has_str:
            reason = "missing_required_cgenff_str_file"
            warnings.append("provide ParamChem/CGenFF .str file; no alternate force field used")
        elif not converter_found:
            reason = "missing_cgenff_charmm2gmx_converter"
            warnings.append("CGenFF stream found but converter is unavailable")
        if not cgenff_found:
            warnings.append("cgenff_executable_not_on_path")
        rows.append(
            {
                "md_candidate_id": row["md_candidate_id"],
                "molecule_id": row["molecule_id"],
                "ligand_smiles": row["standard_smiles"],
                "ligand_input_file": row["ligand_file"],
                "ligand_protonation_state_id": "stage7_prepared_state",
                "ligand_tautomer_state_id": "stage7_prepared_state",
                "cgenff_str_file": str(str_file) if has_str else "",
                "cgenff_version": "unknown" if not cgenff_found else "installed",
                "cgenff_penalty_max": None,
                "cgenff_penalty_mean": None,
                "cgenff_penalty_status": "not_evaluated_missing_str" if not has_str else "pending_parse",
                "charmm2gmx_status": "available" if converter_found else "missing",
                "ligand_itp_file": "",
                "ligand_prm_file": "",
                "ligand_gro_file": "",
                "parameterization_status": status,
                "rejection_reason": reason,
                "warnings_json": json.dumps(warnings),
            }
        )
    out = pd.DataFrame(rows)
    write_table(paths["processed"] / "ligand_parameterization.parquet", out)
    write_table(paths["processed"] / "ligand_parameterization.csv", out)
    return out
