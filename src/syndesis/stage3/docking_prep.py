from __future__ import annotations

import json
import shutil
import subprocess
from pathlib import Path
from typing import Any

import pandas as pd

from syndesis.common.io import ensure_dir, write_table
from syndesis.stage3.schemas import LIGAND_PREP_COLUMNS, RECEPTOR_PREP_COLUMNS


def copy_file(src: str | Path, dst: Path) -> Path:
    ensure_dir(dst.parent)
    shutil.copyfile(src, dst)
    return dst


def prepare_receptors(ensemble: pd.DataFrame, config: dict[str, Any], paths: dict[str, Path]) -> pd.DataFrame:
    rows = []
    engine = config["docking"]["primary_engine"]
    for _, row in ensemble.iterrows():
        rid = row["receptor_id"]
        out = copy_file(row["receptor_file_path"], paths["docking_receptors"] / f"{rid}.pdb")
        docking_format = paths["docking_receptors"] / f"{rid}.pdbqt"
        obabel = config["prep"].get("obabel_path")
        if not obabel:
            raise RuntimeError("Stage 3 receptor PDBQT preparation requires `prep.obabel_path`.")
        ph_target = str(config["prep"].get("ph_target", 7.4))
        completed = subprocess.run(
            [obabel, str(out), "-O", str(docking_format), "-xr", "-p", ph_target, "--partialcharge", "gasteiger"],
            check=False,
            capture_output=True,
            text=True,
            timeout=180,
        )
        if completed.returncode != 0 or not docking_format.exists() or docking_format.stat().st_size == 0:
            raise RuntimeError(f"OpenBabel receptor PDBQT conversion failed for {rid}: {completed.stderr.strip()[:1000]}")
        rows.append({
            "receptor_id": rid,
            "complex_id": row["complex_id"],
            "pdb_id": row["pdb_id"],
            "auth_asym_id": row["auth_asym_id"],
            "receptor_state": row["state_stratum"],
            "input_receptor_file": row["receptor_file_path"],
            "prepared_receptor_file": str(out),
            "docking_format_file": str(docking_format),
            "engine": engine,
            "hydrogens_added_flag": True,
            "protonation_tool": config["prep"]["receptor_prep_tool"],
            "retained_waters_flag": bool(config["prep"].get("retained_waters_default", False)),
            "retained_water_count": int(row.get("pocket_water_count_5a") or 0) if config["prep"].get("retained_waters_default", False) else 0,
            "preparation_warnings_json": json.dumps([]),
            "preparation_status": "prepared_pdbqt",
        })
    frame = pd.DataFrame(rows, columns=RECEPTOR_PREP_COLUMNS)
    write_table(paths["processed"] / "receptor_docking_prep.parquet", frame)
    write_table(paths["processed"] / "receptor_docking_prep.csv", frame)
    return frame


def prepare_ligands(ensemble: pd.DataFrame, config: dict[str, Any], paths: dict[str, Path]) -> pd.DataFrame:
    rows = []
    engine = config["docking"]["primary_engine"]
    for _, row in ensemble.iterrows():
        ligand_id = f"{row['complex_id']}_{row['ligand_comp_id']}".lower()
        immutable = copy_file(row["native_ligand_sdf_path"], paths["native_reference_poses"] / f"{ligand_id}_native.pdb")
        prepared = copy_file(row["native_ligand_sdf_path"], paths["docking_ligands"] / f"{ligand_id}.pdb")
        docking_format = paths["docking_ligands"] / f"{ligand_id}.pdbqt"
        warnings = []
        status = "prepared_copy"
        obabel = config["prep"].get("obabel_path")
        if obabel:
            ph_target = str(config["prep"].get("ph_target", 7.4))
            completed = subprocess.run(
                [obabel, str(prepared), "-O", str(docking_format), "-p", ph_target, "--partialcharge", "gasteiger"],
                check=False,
                capture_output=True,
                text=True,
                timeout=120,
            )
            if completed.returncode == 0 and docking_format.exists():
                status = "prepared_pdbqt"
            else:
                docking_format = prepared
                warnings.append(f"OpenBabel PDBQT conversion failed: {completed.stderr.strip()[:500]}")
        else:
            docking_format = prepared
            warnings.append("OpenBabel path not configured; ligand docking format is copied PDB.")
        rows.append({
            "ligand_id": ligand_id,
            "complex_id": row["complex_id"],
            "pdb_id_native": row["pdb_id"],
            "ligand_comp_id": row["ligand_comp_id"],
            "ligand_instance_id": row["complex_id"].split("_")[-1],
            "native_ligand_file": row["native_ligand_sdf_path"],
            "immutable_reference_pose_file": str(immutable),
            "standard_smiles": None,
            "prepared_smiles": None,
            "protonation_state_id": "native_copy",
            "tautomer_state_id": "native_copy",
            "conformer_id": "native_copy",
            "prepared_ligand_file": str(prepared),
            "docking_format_file": str(docking_format),
            "engine": engine,
            "ligand_prep_tool": config["prep"]["ligand_prep_tool"],
            "charge_model": "gasteiger" if status == "prepared_pdbqt" else config["prep"]["charge_model"],
            "atom_mapping_status": "identity_atom_order_assumed_for_native_copy",
            "native_to_prepared_atom_map_json": json.dumps({}),
            "preparation_warnings_json": json.dumps(warnings),
            "preparation_status": status,
        })
    frame = pd.DataFrame(rows, columns=LIGAND_PREP_COLUMNS)
    write_table(paths["processed"] / "ligand_docking_prep.parquet", frame)
    write_table(paths["processed"] / "ligand_docking_prep.csv", frame)
    return frame
