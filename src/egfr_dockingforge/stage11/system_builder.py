from __future__ import annotations

import json
import os
import shutil
import subprocess
from pathlib import Path
from typing import Any

import pandas as pd
from rdkit import Chem

from egfr_dockingforge.common.io import write_table
from egfr_dockingforge.stage11.mdp_templates import write_mdp_templates


def _run(command: list[str], cwd: Path, log_path: Path, input_text: str | None = None) -> tuple[int, str]:
    log_path.parent.mkdir(parents=True, exist_ok=True)
    env = os.environ.copy()
    result = subprocess.run(command, cwd=cwd, input=input_text, text=True, capture_output=True, env=env)
    text = f"$ {' '.join(command)}\n\nSTDOUT:\n{result.stdout}\n\nSTDERR:\n{result.stderr}\n"
    log_path.write_text(text, encoding="utf-8")
    return result.returncode, text


def _pose_sdf_for_candidate(best_pose_id: str) -> Path:
    # The docked start pose (ProLIF template SDF) lives in the Stage 8 screen for
    # known controls, but in the Stage 9 mini-screen for accepted analogs. Search
    # both so both candidate origins resolve.
    search_dirs = [
        Path("data/processed/stage8/prolif_ligands"),
        Path("data/processed/stage9/stage8_mini_screen/prolif_ligands"),
    ]
    for d in search_dirs:
        path = d / f"{best_pose_id}.pose_template.sdf"
        if path.exists():
            return path
    raise FileNotFoundError(
        f"Missing docked pose SDF for {best_pose_id}; looked in: "
        + ", ".join(str(d) for d in search_dirs)
    )


def _ligand_parameter_dir(output_itp: str) -> Path:
    return Path(output_itp).resolve().parent


def _itp_atom_names(itp_path: Path) -> tuple[str, list[str]]:
    names: list[str] = []
    moleculetype = None
    section = None
    for raw in itp_path.read_text(encoding="utf-8").splitlines():
        line = raw.strip()
        if not line or line.startswith(";"):
            continue
        if line.startswith("[") and line.endswith("]"):
            section = line.strip("[]").strip()
            continue
        if section == "moleculetype" and moleculetype is None:
            moleculetype = line.split()[0]
        elif section == "atoms":
            parts = line.split()
            if len(parts) >= 5 and parts[0].isdigit():
                names.append(parts[4])
    if not moleculetype or not names:
        raise ValueError(f"Could not parse ligand atom names from {itp_path}")
    return moleculetype, names


def _write_pose_ligand_gro(pose_sdf: Path, atom_names: list[str], out_path: Path) -> None:
    mol = Chem.SDMolSupplier(str(pose_sdf), removeHs=False)[0]
    if mol is None:
        raise ValueError(f"RDKit could not read {pose_sdf}")
    mol = Chem.AddHs(mol, addCoords=True)
    if mol.GetNumAtoms() != len(atom_names):
        raise ValueError(f"Pose atom count {mol.GetNumAtoms()} does not match topology atom count {len(atom_names)} for {pose_sdf}")
    conf = mol.GetConformer()
    lines = [f"Docked pose ligand from {pose_sdf.name}", f"{mol.GetNumAtoms():5d}"]
    for idx, atom_name in enumerate(atom_names, start=1):
        pos = conf.GetAtomPosition(idx - 1)
        lines.append(f"{1:5d}{'UNL':>5}{atom_name:>5}{idx:5d}{pos.x/10:8.3f}{pos.y/10:8.3f}{pos.z/10:8.3f}")
    lines.append("   0.00000   0.00000   0.00000")
    out_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _read_gro(path: Path) -> tuple[str, list[str], str]:
    lines = path.read_text(encoding="utf-8").splitlines()
    return lines[0], lines[2:-1], lines[-1]


def _combine_gro(protein_gro: Path, ligand_gro: Path, out_path: Path) -> int:
    _, protein_atoms, box = _read_gro(protein_gro)
    _, ligand_atoms, _ = _read_gro(ligand_gro)
    atoms = protein_atoms + ligand_atoms
    out_path.write_text("Protein-ligand complex\n" + f"{len(atoms):5d}\n" + "\n".join(atoms) + "\n" + box + "\n", encoding="utf-8")
    return len(atoms)


def _insert_ligand_topology(topol_path: Path, ligand_itp: Path, moleculetype: str) -> Path:
    lines = topol_path.read_text(encoding="utf-8").splitlines()
    include_line = f'#include "{ligand_itp.as_posix()}"'
    insert_at = None
    for i, line in enumerate(lines):
        if line.strip().startswith("#include") and "forcefield.itp" in line:
            insert_at = i + 1
            break
    if insert_at is None:
        raise ValueError(f"Could not find forcefield include in {topol_path}")
    lines.insert(insert_at, include_line)
    lines.append(f"{moleculetype:<20} 1")
    topol_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return topol_path


def _water_name(config: dict[str, Any]) -> str:
    water = config["forcefield"]["water_model"].lower()
    if water not in {"opc", "opc3", "tip3p", "spc", "spce", "tip4pew"}:
        raise ValueError(f"Unsupported GROMACS water model for pdb2gmx: {water}")
    return water


def _build_one_system(row: dict[str, Any], param: dict[str, Any], config: dict[str, Any], paths: dict[str, Path]) -> dict[str, Any]:
    gmx = config["forcefield"]["gromacs_executable"]
    work = paths["md_root"] / f"mdsys_{row['md_candidate_id']}"
    if work.exists():
        shutil.rmtree(work)
    work.mkdir(parents=True, exist_ok=True)
    mdp_files = write_mdp_templates(config, work)
    protein_gro = work / "protein.gro"
    protein_top = work / "topol.top"
    posre = work / "posre.itp"
    code, _ = _run(
        [gmx, "pdb2gmx", "-f", row["receptor_file"], "-o", str(protein_gro), "-p", str(protein_top), "-i", str(posre), "-ff", config["forcefield"]["protein_ff"], "-water", _water_name(config), "-ignh"],
        Path.cwd(),
        work / "pdb2gmx.log",
    )
    if code != 0:
        raise RuntimeError(f"pdb2gmx failed for {row['md_candidate_id']}; see {work / 'pdb2gmx.log'}")
    ligand_itp_value = param.get("output_itp") or param.get("ligand_itp_file")
    if not ligand_itp_value:
        raise KeyError("Missing ligand ITP path in parameterization table")
    ligand_itp = _ligand_parameter_dir(ligand_itp_value) / Path(ligand_itp_value).name
    moleculetype, atom_names = _itp_atom_names(ligand_itp)
    ligand_pose_gro = work / "ligand_pose.gro"
    _write_pose_ligand_gro(_pose_sdf_for_candidate(row["best_pose_id"]), atom_names, ligand_pose_gro)
    complex_gro = work / "complex.gro"
    total_atoms = _combine_gro(protein_gro, ligand_pose_gro, complex_gro)
    _insert_ligand_topology(protein_top, ligand_itp, moleculetype)
    boxed = work / "boxed.gro"
    code, _ = _run([gmx, "editconf", "-f", str(complex_gro), "-o", str(boxed), "-bt", config["forcefield"]["box_type"], "-d", str(config["forcefield"]["box_padding_nm"]), "-c"], Path.cwd(), work / "editconf.log")
    if code != 0:
        raise RuntimeError(f"editconf failed for {row['md_candidate_id']}; see {work / 'editconf.log'}")
    solvated = work / "solvated.gro"
    code, _ = _run([gmx, "solvate", "-cp", str(boxed), "-cs", "spc216.gro", "-p", str(protein_top), "-o", str(solvated)], Path.cwd(), work / "solvate.log")
    if code != 0:
        raise RuntimeError(f"solvate failed for {row['md_candidate_id']}; see {work / 'solvate.log'}")
    ions_tpr = work / "ions.tpr"
    code, _ = _run([gmx, "grompp", "-f", str(mdp_files["ions"]), "-c", str(solvated), "-p", str(protein_top), "-o", str(ions_tpr), "-maxwarn", "2"], Path.cwd(), work / "grompp_ions.log")
    if code != 0:
        raise RuntimeError(f"grompp ions failed for {row['md_candidate_id']}; see {work / 'grompp_ions.log'}")
    ionized = work / "ionized.gro"
    code, _ = _run(
        [gmx, "genion", "-s", str(ions_tpr), "-o", str(ionized), "-p", str(protein_top), "-pname", "NA", "-nname", "CL", "-neutral", "-conc", str(config["forcefield"]["salt_concentration_molar"])],
        Path.cwd(),
        work / "genion.log",
        input_text="SOL\n",
    )
    if code != 0:
        raise RuntimeError(f"genion failed for {row['md_candidate_id']}; see {work / 'genion.log'}")
    final_atoms = int(Path(ionized).read_text(encoding="utf-8").splitlines()[1].strip())
    return {
        "md_system_id": f"mdsys_{row['md_candidate_id']}",
        "md_candidate_id": row["md_candidate_id"],
        "protein_file": row["receptor_file"],
        "ligand_file": str(ligand_pose_gro),
        "complex_file": str(complex_gro),
        "topology_file": str(protein_top),
        "position_restraints_file": str(posre),
        "boxed_structure_file": str(boxed),
        "solvated_structure_file": str(solvated),
        "ionized_structure_file": str(ionized),
        "water_model": config["forcefield"]["water_model"],
        "box_type": config["forcefield"]["box_type"],
        "box_padding_nm": float(config["forcefield"]["box_padding_nm"]),
        "num_protein_atoms": int(Path(protein_gro).read_text(encoding="utf-8").splitlines()[1].strip()),
        "num_ligand_atoms": len(atom_names),
        "num_waters": max(final_atoms - total_atoms, 0) // 3,
        "num_na": 0,
        "num_cl": 0,
        "net_charge_before_ions": None,
        "final_charge": None,
        "salt_concentration_molar": float(config["forcefield"]["salt_concentration_molar"]),
        "build_status": "ready_for_md",
        "warnings_json": json.dumps(["solvated_with_spc216_coordinate_box_for_configured_three_site_water_model"]),
    }


def _blocked_system(row: dict, param: dict) -> dict:
    reason = param.get("rejection_reason") or "ligand_parameters_not_ready"
    return {
        "md_system_id": f"mdsys_{row['md_candidate_id']}",
        "md_candidate_id": row["md_candidate_id"],
        "protein_file": row.get("receptor_file", ""),
        "ligand_file": "",
        "complex_file": "",
        "topology_file": "",
        "position_restraints_file": "",
        "boxed_structure_file": "",
        "solvated_structure_file": "",
        "ionized_structure_file": "",
        "water_model": config_get(param, "water_model", ""),
        "box_type": "",
        "box_padding_nm": 0.0,
        "num_protein_atoms": 0,
        "num_ligand_atoms": 0,
        "num_waters": 0,
        "num_na": 0,
        "num_cl": 0,
        "net_charge_before_ions": None,
        "final_charge": None,
        "salt_concentration_molar": 0.0,
        "build_status": "blocked_ligand_parameters_not_ready",
        "warnings_json": json.dumps([f"system_build_blocked:{reason}"]),
    }


def config_get(d: dict, key: str, default):
    return d.get(key, default) if isinstance(d, dict) else default


def build_md_systems(candidates: pd.DataFrame, params: pd.DataFrame, config: dict[str, Any], paths: dict[str, Path]) -> pd.DataFrame:
    pidx = params.set_index("md_candidate_id")
    rows = []
    for row in candidates.to_dict("records"):
        # Only build systems for MD finalists; analogs retained as future work
        # ("selected_md_pending") are skipped here on purpose.
        if "is_md_finalist" in row and not bool(row.get("is_md_finalist", False)):
            continue
        if row["md_candidate_id"] not in pidx.index:
            continue
        param = pidx.loc[row["md_candidate_id"]].to_dict()
        ready = param["parameterization_status"] == "ready"
        if not ready:
            # A finalist whose ligand parameters failed (e.g. charge non-convergence)
            # cannot be simulated honestly; record the block instead of crashing.
            rows.append(_blocked_system(row, param))
            continue
        rows.append(_build_one_system(row, param, config, paths))
    out = pd.DataFrame(rows)
    write_table(paths["processed"] / "system_builds.parquet", out)
    write_table(paths["processed"] / "system_builds.csv", out)
    return out
