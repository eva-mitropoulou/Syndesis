from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pandas as pd

from egfr_dockingforge.common.io import write_table
from egfr_dockingforge.stage5.prolif_engine import has_hydrogens, prepare_ligand_for_prolif, prepare_protein_for_prolif
from egfr_dockingforge.stage5.schemas import INTERACTION_COMPLEX_COLUMNS


def build_interaction_complexes(inputs: dict[str, pd.DataFrame], paths: dict[str, Path]) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    ensemble_ids = set(inputs["ensemble"]["complex_id"].astype(str).str.upper())
    benchmark = inputs["benchmark"]
    for _, row in benchmark[benchmark["include_in_stage1_benchmark"].fillna(False).astype(bool)].iterrows():
        complex_id = str(row["complex_id"])
        if ensemble_ids and complex_id.upper() not in ensemble_ids:
            continue
        source_ligand = row.get("native_ligand_sdf_path")
        ligand = (
            str(prepare_ligand_for_prolif(source_ligand, paths["processed"] / "prolif_ligands"))
            if isinstance(source_ligand, str) and Path(source_ligand).exists()
            else source_ligand
        )
        source_protein = row.get("receptor_clean_path")
        protein = (
            str(prepare_protein_for_prolif(source_protein, paths["processed"] / "prolif_receptors"))
            if isinstance(source_protein, str) and Path(source_protein).exists()
            else source_protein
        )
        rows.append(
            {
                "complex_analysis_id": f"native__{complex_id.lower()}",
                "pose_id": None,
                "complex_id": complex_id,
                "ligand_id": f"{complex_id.lower()}_{str(row.get('ligand_comp_id')).lower()}",
                "receptor_id": complex_id.lower(),
                "pdb_id": row.get("pdb_id"),
                "task_type": "native_cocrystal",
                "is_native_complex": True,
                "protein_file": protein,
                "ligand_file": ligand,
                "complex_file": row.get("native_complex_path"),
                "ligand_source_coordinate_type": "native_ligand",
                "hydrogens_present_flag": has_hydrogens(ligand) if isinstance(ligand, str) and Path(ligand).exists() else False,
                "ligand_bond_orders_available_flag": False,
                "protonation_state_id": "native",
                "tautomer_state_id": "native",
                "complex_build_status": "ready" if isinstance(ligand, str) and Path(ligand).exists() and isinstance(protein, str) and Path(protein).exists() else "missing_input_file",
                "warnings_json": json.dumps([]),
            }
        )
    poses = inputs["docked_poses"]
    scores = inputs["pose_scores"][["pose_id", "task_type"]]
    merged = poses.merge(scores, on="pose_id", how="left")
    for _, row in merged.iterrows():
        source_protein = Path("data/processed/stage3/docking_receptors") / f"{row['target_receptor_id']}.pdb"
        protein = prepare_protein_for_prolif(source_protein, paths["processed"] / "prolif_receptors") if source_protein.exists() else source_protein
        source_ligand = row.get("pose_file")
        ligand = (
            str(prepare_ligand_for_prolif(source_ligand, paths["processed"] / "prolif_ligands"))
            if isinstance(source_ligand, str) and Path(source_ligand).exists()
            else source_ligand
        )
        rows.append(
            {
                "complex_analysis_id": f"pose__{row['pose_id']}",
                "pose_id": row["pose_id"],
                "complex_id": None,
                "ligand_id": row["ligand_id"],
                "receptor_id": row["target_receptor_id"],
                "pdb_id": str(row["target_receptor_id"]).split("_")[0].upper(),
                "task_type": row.get("task_type"),
                "is_native_complex": False,
                "protein_file": str(protein),
                "ligand_file": ligand,
                "complex_file": None,
                "ligand_source_coordinate_type": "stage3_docked_pose",
                "hydrogens_present_flag": has_hydrogens(ligand) if isinstance(ligand, str) and Path(ligand).exists() else False,
                "ligand_bond_orders_available_flag": False,
                "protonation_state_id": row.get("protonation_state_id"),
                "tautomer_state_id": row.get("tautomer_state_id"),
                "complex_build_status": "ready" if isinstance(ligand, str) and Path(ligand).exists() and protein.exists() else "missing_input_file",
                "warnings_json": json.dumps([]),
            }
        )
    frame = pd.DataFrame(rows, columns=INTERACTION_COMPLEX_COLUMNS)
    write_table(paths["processed"] / "interaction_complexes.parquet", frame)
    write_table(paths["processed"] / "interaction_complexes.csv", frame)
    return frame
