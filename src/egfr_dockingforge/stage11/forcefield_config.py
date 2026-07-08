from __future__ import annotations

import json
import shutil
import subprocess
from pathlib import Path
from typing import Any

from egfr_dockingforge.common.io import resolve_path, write_json


def gromacs_version(executable: str) -> str:
    if not Path(executable).exists() and not shutil.which(executable):
        return "missing"
    try:
        completed = subprocess.run([executable, "--version"], capture_output=True, text=True, timeout=20, check=False)
    except OSError:
        return "missing"
    for line in completed.stdout.splitlines():
        if "GROMACS version" in line:
            return line.split(":", 1)[-1].strip()
    return completed.stdout.splitlines()[0].strip() if completed.stdout else "unknown"


def tool_version(executable: str, args: list[str] | None = None) -> str:
    path = resolve_path(executable)
    cmd = [str(path) if path.exists() else executable] + (args or ["--version"])
    try:
        completed = subprocess.run(cmd, capture_output=True, text=True, timeout=30, check=False)
    except OSError:
        return "missing"
    text = (completed.stdout + "\n" + completed.stderr).strip()
    for line in text.splitlines():
        clean = line.strip()
        if "ACPYPE:" in clean:
            return clean
        if "AmberTools" in clean or "antechamber" in clean.lower():
            return clean
    return text.splitlines()[0].strip() if text else "unknown"


def ambertools_version_from_prefix(prefix: str) -> str:
    path = resolve_path(prefix) / "conda-meta"
    if path.exists():
        matches = sorted(path.glob("ambertools-*.json"))
        if matches:
            name = matches[0].name.removesuffix(".json")
            return name
    return "unknown"


def acpype_version_from_prefix(prefix: str) -> str:
    path = resolve_path(prefix) / "conda-meta"
    if path.exists():
        matches = sorted(path.glob("acpype-*.json"))
        if matches:
            name = matches[0].name.removesuffix(".json")
            return name
    return "unknown"


def validate_forcefield_stack(config: dict[str, Any]) -> list[str]:
    ff = config["forcefield"]
    warnings = []
    backend = ff.get("ligand_parameterization_backend", "amber_gaff2_acpype")
    if backend == "amber_gaff2_acpype":
        if not ff["protein_ff"].lower().startswith("amber"):
            warnings.append("gaff_ligand_with_non_amber_protein_forcefield")
        if ff["ligand_ff"].lower() != "gaff2":
            warnings.append("non_gaff2_ligand_forcefield_for_amber_backend")
    elif backend == "charmm_cgenff_paramchem":
        if ff["protein_ff"] != "CHARMM36m" or ff["ligand_ff"] != "CGenFF" or "TIP3P" not in ff["water_model"]:
            warnings.append("incompatible_or_nondefault_charmm_cgenff_stack")
    else:
        warnings.append(f"unknown_ligand_parameterization_backend:{backend}")
    return warnings


def write_forcefield_policy(config: dict[str, Any], paths: dict[str, Path]) -> dict:
    ff = config["forcefield"]
    protein_ff_path = ""
    if ff["protein_ff"].lower().startswith("amber"):
        candidate = Path("/usr/local/gromacs/share/gromacs/top") / f"{ff['protein_ff'].lower()}.ff"
        protein_ff_path = str(candidate) if candidate.exists() else ""
    else:
        protein_ff_path = ff.get("charmm36_path", "")
    policy = {
        "ligand_parameterization_backend": ff.get("ligand_parameterization_backend", "amber_gaff2_acpype"),
        "protein_force_field": ff["protein_ff"],
        "protein_force_field_path": protein_ff_path,
        "ligand_force_field_method": ff["ligand_ff"],
        "ligand_force_field_version": "GAFF2 via AmberTools/ACPYPE" if ff.get("ligand_parameterization_backend") == "amber_gaff2_acpype" else "configured_CGenFF",
        "ligand_charge_model": ff.get("ligand_charge_model", ""),
        "water_model": ff["water_model"],
        "ion_parameters": ff["ion_model"],
        "gromacs_executable": ff["gromacs_executable"],
        "gromacs_version": gromacs_version(ff["gromacs_executable"]),
        "ambertools_version": ambertools_version_from_prefix(ff.get("ambertools_prefix", "")),
        "acpype_version": acpype_version_from_prefix(ff.get("ambertools_prefix", "")),
        "ligand_charge_protonation_source": "previous_stage_prepared_ligand_state",
        "parameterization_warnings": validate_forcefield_stack(config),
        "exact_files_used": {},
    }
    write_json(paths["processed"] / "forcefield_policy.json", policy)
    return policy
