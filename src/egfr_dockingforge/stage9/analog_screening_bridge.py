from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pandas as pd
import yaml
from rdkit import Chem
from rdkit.Chem import AllChem

from egfr_dockingforge.common.io import ensure_dir, load_yaml, project_root, resolve_path, write_table
from egfr_dockingforge.stage8.candidate_screening import run_stage8_all


def _write_sdf(smiles: str, analog_id: str, out_dir: Path) -> str:
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    if mol is None:
        raise ValueError(f"Invalid analog SMILES for SDF: {analog_id}")
    AllChem.EmbedMolecule(mol, randomSeed=909)
    AllChem.UFFOptimizeMolecule(mol, maxIters=200)
    target = out_dir / f"{analog_id}.sdf"
    writer = Chem.SDWriter(str(target))
    mol.SetProp("_Name", analog_id)
    writer.write(mol)
    writer.close()
    return str(target)


def write_stage8_mini_input(
    candidates: pd.DataFrame,
    validation: pd.DataFrame,
    config: dict[str, Any],
    paths: dict[str, Path],
) -> tuple[pd.DataFrame, pd.DataFrame, Path]:
    valid_ids = set(validation[validation["hard_scope_pass"].astype(bool)]["analog_id"])
    selected = candidates[candidates["analog_id"].isin(valid_ids)].head(int(config["screening"]["max_analogs_to_screen"])).copy()
    rows = []
    master_rows = []
    for row in selected.to_dict("records"):
        sdf = _write_sdf(row["standard_smiles"], row["analog_id"], paths["analog_sdf_dir"])
        prep_id = f"prep_{row['analog_id']}"
        rows.append(
            {
                "prepared_ligand_id": prep_id,
                "molecule_id": row["analog_id"],
                "source": "stage9_analog",
                "screening_role": "stage9_analog",
                "screening_subset": row["strategy_name"],
                "standard_smiles": row["standard_smiles"],
                "prepared_smiles": row["standard_smiles"],
                "sdf_path": sdf,
                "novelty_bucket": "stage9_analog",
                "closest_known_molecule_id": row["parent_molecule_id"],
                "tanimoto_to_closest_known": row["parent_tanimoto"],
            }
        )
        master_rows.append(
            {
                "molecule_id": row["analog_id"],
                "subsource": row["strategy_name"],
                "hard_scope_pass": True,
                "include_in_screening_library": True,
                "risk_flags_json": "[]",
                "scaffold_id": row["seed_id"],
            }
        )
    stage8_input = pd.DataFrame(rows)
    master = pd.DataFrame(master_rows)
    write_table(paths["processed"] / "stage8_analog_screening_input.parquet", stage8_input)
    write_table(paths["processed"] / "stage8_analog_master.parquet", master)
    stage8_cfg = load_yaml(resolve_path(config["inputs"]["stage8_config"]))
    stage8_cfg["inputs"]["stage7_screening_input"] = str(paths["processed"] / "stage8_analog_screening_input.parquet")
    stage8_cfg["inputs"]["candidate_library_master"] = str(paths["processed"] / "stage8_analog_master.parquet")
    stage8_cfg["paths"]["processed"] = config["screening"]["stage8_processed"]
    stage8_cfg["screening"]["max_known_controls"] = 0
    stage8_cfg["screening"]["max_candidates"] = int(config["screening"]["max_analogs_to_screen"])
    stage8_cfg["screening"]["receptor_limit"] = int(config["screening"]["receptor_limit"])
    stage8_cfg["docking"]["num_modes"] = int(config["screening"]["num_modes"])
    cfg_path = resolve_path(config["screening"]["stage8_mini_config"], project_root())
    ensure_dir(cfg_path.parent)
    cfg_path.write_text(yaml.safe_dump(stage8_cfg, sort_keys=False), encoding="utf-8")
    return stage8_input, master, cfg_path


def screen_analog_batch(candidates: pd.DataFrame, validation: pd.DataFrame, seeds: pd.DataFrame, config: dict[str, Any], paths: dict[str, Path]) -> pd.DataFrame:
    if not bool(config["screening"]["enabled"]):
        raise RuntimeError("Stage 9 screening is disabled; prompt requires Stage 8 mini-screen.")
    _, _, mini_config = write_stage8_mini_input(candidates, validation, config, paths)
    run_stage8_all(mini_config)
    stage8_processed = resolve_path(config["screening"]["stage8_processed"], project_root())
    agg = pd.read_parquet(stage8_processed / "candidate_aggregate_scores.parquet")
    candidate_idx = candidates.set_index("analog_id")
    seed_idx = seeds.set_index("seed_id")
    rows = []
    for row in agg.to_dict("records"):
        analog = candidate_idx.loc[row["molecule_id"]]
        seed = seed_idx.loc[analog["seed_id"]]
        heavy = Chem.MolFromSmiles(analog["standard_smiles"]).GetNumHeavyAtoms()
        le = -float(row["best_gnina_cnnaffinity"]) / max(heavy, 1) if pd.notna(row["best_gnina_cnnaffinity"]) else 0.0
        binding_preserved = (
            float(row["best_pose_confidence"]) >= float(config["acceptance"]["min_pose_confidence"])
            and float(row["best_key_interaction_recall_consensus"]) >= float(config["acceptance"]["min_key_interaction_recall"])
        )
        rows.append(
            {
                "analog_id": row["molecule_id"],
                "seed_id": analog["seed_id"],
                "strategy_name": analog["strategy_name"],
                "iteration_id": analog["iteration_id"],
                "prepared_ligand_id": f"prep_{row['molecule_id']}",
                "best_pose_id": row["best_screening_pose_id"],
                "best_receptor_id": row["best_target_receptor_id"],
                "best_receptor_state": row["best_receptor_state"],
                "best_docking_score": row["best_docking_score"],
                "best_gnina_cnnscore": row["best_gnina_cnnscore"],
                "best_gnina_cnnaffinity": row["best_gnina_cnnaffinity"],
                "best_ifp_tanimoto_to_consensus": row["best_ifp_tanimoto_to_consensus"],
                "best_key_interaction_recall_consensus": row["best_key_interaction_recall_consensus"],
                "best_pose_confidence": row["best_pose_confidence"],
                "best_calibrated_confidence": row["best_calibrated_confidence"],
                "ligand_efficiency": le,
                "binding_mode_preserved_flag": bool(binding_preserved),
                "stage8_screening_status": row["candidate_decision_label"],
                "warnings_json": json.dumps(["parent_seed_rank_%s" % seed["parent_candidate_rank"]]),
            }
        )
    out = pd.DataFrame(rows)
    write_table(paths["processed"] / "analog_screening_results.parquet", out)
    write_table(paths["processed"] / "analog_screening_results.csv", out)
    return out
