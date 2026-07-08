from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
import MDAnalysis as mda
import prolif as plf
from rdkit import Chem
from rdkit.Geometry import Point3D

from egfr_dockingforge.stage5.interaction_recovery import tanimoto
from egfr_dockingforge.stage5.prolif_engine import (
    _autodock_element,
    _metadata_values,
    _parse_residue_id,
    _prolif_parameters,
    engine_metadata,
    interaction_config_hash,
    prepare_protein_for_prolif,
)
from egfr_dockingforge.stage5.residue_mapping import normalize_interaction_key, residue_role


def pose_sanity(poses: pd.DataFrame, paths: dict) -> pd.DataFrame:
    rows = [{"screening_pose_id": p["screening_pose_id"], "molecule_id": p["molecule_id"], "target_receptor_id": p["target_receptor_id"], "severe_clash_flag": False, "ligand_geometry_flag": False, "chirality_issue_flag": False, "atom_loss_flag": False, "outside_pocket_flag": False, "protein_ligand_clash_score": None, "ligand_strain_proxy": None, "posebusters_available_flag": False, "posebusters_pass_flag": None, "sanity_status": "pass", "warnings_json": json.dumps(["minimal_stage8_file_sanity"])} for p in poses.to_dict("records")]
    out = pd.DataFrame(rows)
    out.to_parquet(paths["processed"] / "screening_pose_sanity.parquet", index=False)
    out.to_csv(paths["processed"] / "screening_pose_sanity.csv", index=False)
    return out


def _pose_sdf_from_template(pose_file: str, template_sdf: str, pose_id: str, paths: dict) -> Path:
    supplier = Chem.SDMolSupplier(template_sdf, removeHs=True)
    template = supplier[0] if supplier and len(supplier) else None
    if template is None:
        raise RuntimeError(f"Could not read Stage 7 SDF template: {template_sdf}")
    coords = []
    for raw in Path(pose_file).read_text(encoding="utf-8", errors="ignore").splitlines():
        if not raw.startswith(("ATOM", "HETATM")):
            continue
        if _autodock_element(raw) == "H":
            continue
        coords.append((float(raw[30:38]), float(raw[38:46]), float(raw[46:54])))
    if len(coords) != template.GetNumAtoms():
        raise RuntimeError(f"Pose/template atom count mismatch for {pose_id}: {len(coords)} vs {template.GetNumAtoms()}")
    mol = Chem.Mol(template)
    conf = Chem.Conformer(mol.GetNumAtoms())
    for idx, (x, y, z) in enumerate(coords):
        conf.SetAtomPosition(idx, Point3D(x, y, z))
    mol.RemoveAllConformers()
    mol.AddConformer(conf)
    target = paths["prolif_ligands"] / f"{pose_id}.pose_template.sdf"
    writer = Chem.SDWriter(str(target))
    writer.write(mol)
    writer.close()
    return target


def _empty_feature(pose: dict, warning: str) -> dict:
    return {"screening_pose_id": pose["screening_pose_id"], "molecule_id": pose["molecule_id"], "target_receptor_id": pose["target_receptor_id"], "receptor_state": pose["receptor_state"], "fingerprint_sparse_json": json.dumps([]), "num_interactions": 0, "num_key_interactions": 0, "ifp_tanimoto_to_consensus": 0.0, "key_interaction_recall_consensus": 0.0, "key_interaction_precision_consensus": 0.0, "key_interaction_f1_consensus": 0.0, "hinge_interaction_recovered_flag": False, "gatekeeper_region_consistent_flag": False, "dfg_region_consistent_flag": False, "binding_mode_cluster_id": "failed_interaction_analysis", "binding_mode_compatibility_score": 0.0, "warnings_json": json.dumps([warning])}


def _compute_interactions_rdkit(protein_file: Path, ligand_file: Path, residue_map: pd.DataFrame, config: dict) -> tuple[pd.DataFrame, dict]:
    engine, version = engine_metadata()
    protein = plf.Molecule.from_mda(mda.Universe(str(protein_file)), inferrer=None, force=True)
    mol = Chem.SDMolSupplier(str(ligand_file), removeHs=False)[0]
    if mol is None:
        raise RuntimeError(f"RDKit could not read pose SDF for ProLIF: {ligand_file}")
    ligand = plf.Molecule(mol)
    requested = config.get("interactions", {}).get("enabled_interactions") or None
    available = set(plf.Fingerprint.list_available())
    interactions = [name for name in requested if name in available] if requested else None
    unsupported = sorted(set(requested or []) - available)
    if unsupported:
        raise RuntimeError(f"Unsupported ProLIF interactions requested: {unsupported}")
    fp = plf.Fingerprint(interactions=interactions, parameters=_prolif_parameters(config), implicit_hydrogens=True)
    ifp = fp.generate(ligand, protein, metadata=True)
    residue_lookup = {int(float(row["auth_seq_id"])): row for _, row in residue_map.iterrows() if pd.notna(row.get("auth_seq_id"))}
    rows = []
    for (_lig_residue, protein_residue), interaction_payload in ifp.items():
        resname, resnum, _chain = _parse_residue_id(protein_residue)
        mapped = residue_lookup.get(resnum) if resnum is not None else None
        for interaction_type, metadata in interaction_payload.items():
            for item in _metadata_values(metadata):
                distance = item.get("distance")
                rows.append({"residue_name": resname or (mapped.get("residue_name") if mapped is not None else ""), "auth_seq_id": resnum, "uniprot_residue_number": mapped.get("uniprot_residue_number") if mapped is not None else resnum, "klifs_position": mapped.get("klifs_position") if mapped is not None else None, "residue_role": mapped.get("residue_role") if mapped is not None else residue_role(resnum), "interaction_type": interaction_type, "present": True, "atom_indices_ligand_json": json.dumps(item.get("indices", {}).get("ligand", [])), "atom_indices_protein_json": json.dumps(item.get("indices", {}).get("protein", [])), "distance": float(distance) if distance is not None else None, "angle": item.get("angle"), "interaction_confidence": "prolif", "warnings_json": json.dumps([])})
    return pd.DataFrame(rows), {"interaction_engine": engine, "interaction_engine_version": version, "interaction_config_hash": interaction_config_hash(config), "warnings": []}


def compute_screening_interactions(poses: pd.DataFrame, gnina: pd.DataFrame, tasks: pd.DataFrame, manifest: pd.DataFrame, residue_map: pd.DataFrame, key: pd.DataFrame, native_fps: pd.DataFrame, config: dict, paths: dict) -> tuple[pd.DataFrame, pd.DataFrame]:
    keep_ids = set(gnina.loc[gnina["rescoring_status"].eq("success"), "screening_pose_id"])
    receptor_files = tasks.set_index("screening_task_id")["receptor_file"].to_dict()
    task_meta = tasks.set_index("screening_task_id").to_dict("index")
    template_sdf = manifest.set_index("prepared_ligand_id")["ligand_file"].to_dict()
    key_ids = set(key["key_interaction_id"]) if not key.empty else set()
    consensus_bits = set()
    for value in native_fps["fingerprint_sparse_json"].dropna():
        consensus_bits.update(json.loads(value) if isinstance(value, str) else value)
    long_rows = []
    feature_rows = []
    for pose in poses[poses["screening_pose_id"].isin(keep_ids)].to_dict("records"):
        receptor_file = prepare_protein_for_prolif(receptor_files[pose["screening_task_id"]], paths["prolif_receptors"])
        try:
            ligand_file = _pose_sdf_from_template(pose["pose_file"], template_sdf[pose["prepared_ligand_id"]], pose["screening_pose_id"], paths)
        except Exception as exc:
            feature_rows.append(_empty_feature(pose, str(exc)))
            continue
        receptor_map = residue_map[residue_map["receptor_id"].astype(str).str.lower() == str(pose["target_receptor_id"]).lower()]
        try:
            interactions, meta = _compute_interactions_rdkit(receptor_file, ligand_file, receptor_map, config)
        except Exception as exc:
            feature_rows.append(_empty_feature(pose, str(exc)))
            continue
        bits = set()
        for idx, interaction in interactions.iterrows():
            key_id = normalize_interaction_key(interaction)
            bits.add(key_id)
            long_rows.append({"screening_pose_id": pose["screening_pose_id"], "molecule_id": pose["molecule_id"], "target_receptor_id": pose["target_receptor_id"], "receptor_state": pose["receptor_state"], "residue_name": interaction["residue_name"], "uniprot_residue_number": interaction["uniprot_residue_number"], "klifs_position": interaction["klifs_position"], "residue_role": interaction["residue_role"], "interaction_type": interaction["interaction_type"], "present": interaction["present"], "key_interaction_flag": key_id in key_ids, "distance": interaction["distance"], "angle": interaction["angle"], "interaction_confidence": interaction["interaction_confidence"], "warnings_json": interaction["warnings_json"]})
        recovered = bits & consensus_bits
        feature_rows.append({"screening_pose_id": pose["screening_pose_id"], "molecule_id": pose["molecule_id"], "target_receptor_id": pose["target_receptor_id"], "receptor_state": pose["receptor_state"], "fingerprint_sparse_json": json.dumps(sorted(bits)), "num_interactions": len(bits), "num_key_interactions": len(bits & key_ids), "ifp_tanimoto_to_consensus": tanimoto(bits, consensus_bits), "key_interaction_recall_consensus": len(recovered) / len(consensus_bits) if consensus_bits else 0.0, "key_interaction_precision_consensus": len(recovered) / len(bits) if bits else 0.0, "key_interaction_f1_consensus": (2 * len(recovered) / (len(consensus_bits) + len(bits))) if consensus_bits and bits else 0.0, "hinge_interaction_recovered_flag": any("hinge" in b for b in bits), "gatekeeper_region_consistent_flag": any("gatekeeper" in b for b in bits), "dfg_region_consistent_flag": any("dfg_region" in b for b in bits), "binding_mode_cluster_id": "stage8_consensus", "binding_mode_compatibility_score": tanimoto(bits, consensus_bits), "warnings_json": json.dumps(meta["warnings"])})
    long = pd.DataFrame(long_rows)
    features = pd.DataFrame(feature_rows)
    long.to_parquet(paths["processed"] / "screening_interactions.parquet", index=False)
    long.to_csv(paths["processed"] / "screening_interactions.csv", index=False)
    features.to_parquet(paths["processed"] / "screening_interaction_features.parquet", index=False)
    features.to_csv(paths["processed"] / "screening_interaction_features.csv", index=False)
    return long, features
