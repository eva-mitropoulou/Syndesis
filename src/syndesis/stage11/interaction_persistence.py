from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

from syndesis.common.io import write_table


def _json_list(value) -> list:
    if value is None or (isinstance(value, float) and pd.isna(value)):
        return []
    if isinstance(value, list):
        return value
    try:
        parsed = json.loads(value)
        return parsed if isinstance(parsed, list) else [str(parsed)]
    except Exception:
        return [str(value)]


def _interaction_cutoff(interaction_type: str) -> float:
    lowered = str(interaction_type).lower()
    if "hba" in lowered or "hbd" in lowered or "hydrogen" in lowered:
        return 3.5
    if "salt" in lowered or "ionic" in lowered:
        return 4.0
    return 4.5


def _candidate_residue_ids(row: dict) -> list[int]:
    residue_number = row.get("uniprot_residue_number")
    if residue_number is None or pd.isna(residue_number):
        return []
    uniprot = int(residue_number)
    ids = [uniprot]
    klifs = str(row.get("klifs_position") or "")
    if "offset_-24" in klifs:
        ids.insert(0, uniprot - 24)
    if uniprot - 24 not in ids:
        ids.append(uniprot - 24)
    return ids


def _select_residue(universe, row: dict):
    expected = str(row.get("residue_name") or "").upper()
    protein = universe.select_atoms("protein")
    by_id = {int(res.resid): res for res in protein.residues}
    candidates = [by_id[resid] for resid in _candidate_residue_ids(row) if resid in by_id]
    for residue in candidates:
        if not expected or str(residue.resname).upper().endswith(expected):
            return residue, []
    if candidates:
        residue = candidates[0]
        return residue, [f"residue_name_mismatch_used_resid_{int(residue.resid)}_{residue.resname}_expected_{expected}"]
    if expected:
        uniprot = int(row.get("uniprot_residue_number"))
        nearby = [
            res
            for res in protein.residues
            if abs(int(res.resid) - uniprot) <= 30 and str(res.resname).upper().endswith(expected)
        ]
        if nearby:
            residue = nearby[0]
            return residue, [f"residue_id_offset_inferred_used_resid_{int(residue.resid)}"]
    return None, ["interaction_residue_not_found_in_md_topology"]


def _analyze_interactions_for_metric(metric: dict, keys: list[dict], gmx: str) -> tuple[list[dict], dict]:
    import MDAnalysis as mda
    from MDAnalysis.lib.distances import distance_array

    from syndesis.stage11.trajectory_analysis import pbc_corrected_trajectory

    trajectory = Path(metric["trajectory_file"])
    if not trajectory.exists():
        raise FileNotFoundError("production trajectory missing")
    tpr = trajectory.with_suffix(".tpr")
    if not tpr.exists():
        raise FileNotFoundError("production TPR missing (needed for PBC correction)")
    # Reuse the PBC-corrected, backbone-fitted, water-stripped trajectory produced
    # by the trajectory analysis so interaction distances are measured in the
    # protein frame with periodic images resolved.
    fitted, _ = pbc_corrected_trajectory(gmx, tpr, trajectory, trajectory.parent)
    # Water-stripped topology matches the fitted trajectory's atom set.
    topology = trajectory.parent / (trajectory.stem + "_noPBCwater.gro")
    if not topology.exists():
        topology = trajectory.with_suffix(".gro")
    if not topology.exists():
        topology = trajectory.parent / "production_quick.gro"
    universe = mda.Universe(str(topology), str(fitted))
    ligand = universe.select_atoms("resname UNL and not name H*")
    if len(ligand) == 0:
        ligand = universe.select_atoms("not protein and not resname SOL NA CL and not name H*")
    if len(ligand) == 0:
        raise ValueError("ligand atom selection is empty")

    prepared = []
    for i, key in enumerate(keys, start=1):
        residue, warnings = _select_residue(universe, key)
        cutoff = _interaction_cutoff(str(key.get("interaction_type", "")))
        if residue is None:
            prepared.append((key, i, None, cutoff, warnings))
            continue
        atoms = residue.atoms.select_atoms("not name H*")
        if len(atoms) == 0:
            prepared.append((key, i, None, cutoff, warnings + ["interaction_residue_has_no_heavy_atoms"]))
        else:
            prepared.append((key, i, atoms, cutoff, warnings))

    distances: list[list[float]] = [[] for _ in prepared]
    present: list[list[bool]] = [[] for _ in prepared]
    for ts in universe.trajectory:
        ligand_positions = ligand.positions
        box = ts.dimensions
        for idx, (_, _, atoms, cutoff, _) in enumerate(prepared):
            if atoms is None:
                continue
            # Minimum-image convention via the box so periodic images do not
            # inflate residue-ligand contact distances.
            d = distance_array(atoms.positions, ligand_positions, box=box)
            minimum = float(np.min(d))
            distances[idx].append(minimum)
            present[idx].append(minimum <= cutoff)

    rows = []
    for idx, (key, i, atoms, cutoff, warnings) in enumerate(prepared):
        values = distances[idx]
        hits = present[idx]
        computed = bool(values)
        rows.append(
            {
                "md_candidate_id": metric["md_candidate_id"],
                "md_system_id": metric["md_system_id"],
                "replicate_id": metric["replicate_id"],
                "key_interaction_id": key.get("key_interaction_id", f"key_{i}"),
                "residue_name": key.get("residue_name", ""),
                "uniprot_residue_number": key.get("uniprot_residue_number", ""),
                "klifs_position": key.get("klifs_position", ""),
                "residue_role": key.get("residue_role", ""),
                "interaction_type": key.get("interaction_type", ""),
                "initial_pose_present_flag": bool(hits[0]) if hits else False,
                "consensus_key_interaction_flag": True,
                "occupancy_fraction": float(np.mean(hits)) if hits else 0.0,
                "mean_distance": float(np.mean(values)) if values else None,
                "median_distance": float(np.median(values)) if values else None,
                "persistence_status": "computed" if computed else "not_computed_missing_residue_or_ligand",
                "warnings_json": json.dumps(warnings),
            }
        )
    long = pd.DataFrame(rows)
    valid = long[long["persistence_status"].eq("computed")]
    if valid.empty:
        summary = {
            "ifp_tanimoto_to_initial_median": 0.0,
            "ifp_tanimoto_to_consensus_median": 0.0,
            "key_interaction_occupancy_mean": 0.0,
            "key_interaction_occupancy_min": 0.0,
            "hinge_interaction_occupancy": 0.0,
            "gatekeeper_region_persistence": 0.0,
            "dfg_region_persistence": 0.0,
            "binding_mode_persistence_score": 0.0,
            "binding_mode_preserved_flag": False,
            "warnings_json": json.dumps(["no_interaction_persistence_values_computed"]),
        }
    else:
        role = valid.groupby("residue_role")["occupancy_fraction"].mean().to_dict()
        occupancy_mean = float(valid["occupancy_fraction"].mean())
        occupancy_min = float(valid["occupancy_fraction"].min())
        summary = {
            "ifp_tanimoto_to_initial_median": occupancy_mean,
            "ifp_tanimoto_to_consensus_median": occupancy_mean,
            "key_interaction_occupancy_mean": occupancy_mean,
            "key_interaction_occupancy_min": occupancy_min,
            "hinge_interaction_occupancy": float(role.get("hinge", 0.0)),
            "gatekeeper_region_persistence": float(role.get("gatekeeper", 0.0)),
            "dfg_region_persistence": float(role.get("dfg_region", 0.0)),
            "binding_mode_persistence_score": occupancy_mean,
            "binding_mode_preserved_flag": bool(occupancy_mean >= 0.5 and occupancy_min >= 0.2),
            "warnings_json": json.dumps([]),
        }
    return rows, summary


def compute_interaction_persistence(candidates: pd.DataFrame, metrics: pd.DataFrame, key_interactions: pd.DataFrame | None, paths: dict[str, Path], config: dict) -> tuple[pd.DataFrame, pd.DataFrame]:
    gmx = config["forcefield"]["gromacs_executable"]
    rows = []
    summary_rows = []
    # Persistence is scored against the CONSERVED CORE of key interactions
    # (hinge, gatekeeper, catalytic Lys/Glu, DFG), not the union of every
    # native's idiosyncratic contacts, so that occupancy_min is not dragged to
    # zero by contacts that no single ligand is expected to make. Mirrors the
    # Stage 5 interaction-recovery definition.
    core_min_freq = float(config.get("stability", {}).get("consensus_core_min_frequency", 0.6))
    if key_interactions is not None and not key_interactions.empty and "native_frequency" in key_interactions.columns:
        core = key_interactions[
            (key_interactions["native_frequency"].astype(float) >= core_min_freq)
            | key_interactions.get("manual_override_flag", pd.Series(False, index=key_interactions.index)).fillna(False).astype(bool)
        ]
        keys = (core if not core.empty else key_interactions).to_dict("records")
    else:
        keys = key_interactions.to_dict("records") if key_interactions is not None and not key_interactions.empty else []
    for metric in metrics.to_dict("records"):
        complete = metric.get("trajectory_analysis_status") == "complete"
        if complete:
            try:
                metric_rows, summary = _analyze_interactions_for_metric(metric, keys, gmx)
                rows.extend(metric_rows)
                summary_rows.append(
                    {
                        "md_candidate_id": metric["md_candidate_id"],
                        "md_system_id": metric["md_system_id"],
                        "replicate_id": metric["replicate_id"],
                        **summary,
                    }
                )
                continue
            except Exception as exc:
                metric_warning = f"interaction_persistence_failed:{exc}"
        else:
            metric_warning = "production_trajectory_missing"
        for i, key in enumerate(keys, start=1):
            rows.append(
                {
                    "md_candidate_id": metric["md_candidate_id"],
                    "md_system_id": metric["md_system_id"],
                    "replicate_id": metric["replicate_id"],
                    "key_interaction_id": key.get("key_interaction_id", f"key_{i}"),
                    "residue_name": key.get("residue_name", ""),
                    "uniprot_residue_number": key.get("uniprot_residue_number", ""),
                    "klifs_position": key.get("klifs_position", ""),
                    "residue_role": key.get("residue_role", ""),
                    "interaction_type": key.get("interaction_type", ""),
                    "initial_pose_present_flag": False,
                    "consensus_key_interaction_flag": True,
                    "occupancy_fraction": 0.0,
                    "mean_distance": None,
                    "median_distance": None,
                    "persistence_status": "not_computed_missing_trajectory" if not complete else "failed_analysis_error",
                    "warnings_json": json.dumps([metric_warning]),
                }
            )
        summary_rows.append(
            {
                "md_candidate_id": metric["md_candidate_id"],
                "md_system_id": metric["md_system_id"],
                "replicate_id": metric["replicate_id"],
                "ifp_tanimoto_to_initial_median": 0.0,
                "ifp_tanimoto_to_consensus_median": 0.0,
                "key_interaction_occupancy_mean": 0.0,
                "key_interaction_occupancy_min": 0.0,
                "hinge_interaction_occupancy": 0.0,
                "gatekeeper_region_persistence": 0.0,
                "dfg_region_persistence": 0.0,
                "binding_mode_persistence_score": 0.0,
                "binding_mode_preserved_flag": False,
                "warnings_json": json.dumps([metric_warning]),
            }
        )
    long = pd.DataFrame(rows)
    summary = pd.DataFrame(summary_rows)
    write_table(paths["processed"] / "md_interaction_persistence.parquet", long)
    write_table(paths["processed"] / "md_interaction_persistence.csv", long)
    write_table(paths["processed"] / "md_binding_mode_persistence.parquet", summary)
    write_table(paths["processed"] / "md_binding_mode_persistence.csv", summary)
    return long, summary
