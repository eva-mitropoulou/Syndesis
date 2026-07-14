"""Driver: build the LIT-PCBA EGFR enrichment benchmark end to end.

Steps (each idempotent / resumable):
  prepare  : parse LIT-PCBA actives/decoys -> prepare PDBQT (parallel)
  dock     : GPU-batched dock + GNINA rescore across the ensemble (per-receptor ckpt)
  score    : assemble the three matched scoring arms + compute enrichment (EF/BEDROC/
             AUROC with bootstrap CIs), per receptor and ensemble-consensus.

Config-light on purpose: paths are resolved from the project layout + a small tools
dict, so the campaign can be launched by a resilient shell driver like the MD one.
"""
from __future__ import annotations

import json
import shutil
from pathlib import Path
from typing import Any

import pandas as pd

from egfr_dockingforge.common.io import project_root
from egfr_dockingforge.enrichment.ligand_prep import parse_litpcba_smi, prepare_ligands
from egfr_dockingforge.enrichment.orchestrator import run_campaign
from egfr_dockingforge.enrichment.metrics import enrichment_report


# ---- tools (resolved once; the WSL-native docker CLI is required for GNINA) -------
def default_tools(root: Path) -> dict[str, str]:
    required = {name: shutil.which(name) for name in ("unidock", "obabel", "docker")}
    missing = [name for name, path in required.items() if path is None]
    if missing:
        raise RuntimeError(f"Required enrichment executables are unavailable: {', '.join(missing)}")
    return {
        "unidock": required["unidock"],
        "obabel": required["obabel"],
        "docker": required["docker"],
        "gnina_image": "gnina/gnina:latest",
        "docking_receptors_dir": str(root / "data/processed/stage3/docking_receptors"),
    }


def find_litpcba_egfr(extracted_root: Path) -> dict[str, Path]:
    """Locate EGFR actives/inactives SMILES files inside an extracted LIT-PCBA tree.
    LIT-PCBA stores per-target dirs each with actives.smi / inactives.smi."""
    hits: dict[str, Path] = {}
    for p in extracted_root.rglob("*"):
        if not p.is_file():
            continue
        low = p.name.lower()
        parent = p.parent.name.lower()
        if "egfr" not in str(p).lower():
            continue
        if low in {"actives.smi", "active.smi", "actives_final.smi"}:
            hits["actives"] = p
        elif low in {"inactives.smi", "inactive.smi", "decoys.smi", "inactives_final.smi"}:
            hits["decoys"] = p
    return hits


def step_prepare(actives_smi: Path, decoys_smi: Path, prepared_dir: Path, obabel: str, max_workers: int = 12) -> pd.DataFrame:
    records = parse_litpcba_smi(actives_smi, 1) + parse_litpcba_smi(decoys_smi, 0)
    prepared = prepare_ligands(records, prepared_dir, obabel, max_workers=max_workers)
    prepared.to_parquet(prepared_dir.parent / "ligand_prep.parquet", index=False)
    ok = prepared[prepared["prep_status"].isin(["prepared", "cached"])]
    return ok


def _receptor_table_with_pdbqt(root: Path) -> pd.DataFrame:
    ens = pd.read_parquet(root / "data/processed/stage2/receptor_ensemble_v1.parquet")
    prep = pd.read_parquet(root / "data/processed/stage3/receptor_docking_prep.parquet")
    merged = ens.merge(prep[["receptor_id", "docking_format_file"]], on="receptor_id", how="left")
    return merged


def _rank_desc(s: pd.Series) -> pd.Series:
    return (-s).rank(method="average")


def build_arms(master: pd.DataFrame) -> pd.DataFrame:
    """Reduce the per-(ligand,receptor) master to one row per ligand via
    ensemble-consensus (best score across receptors), then define three matched
    scoring arms as ranking signals (higher = more active-like):

      arm_dock  = -best_docking_score           (more negative Vina => better)
      arm_inter = max_pose CNNscore * (1 + native-union recall)
      arm_pose  = pose-confidence probability    (stage6 model; if present)
      arm_gnina = gnina cnnscore                 (reference CNN arm)

    The interaction arm requires neural scores and strict same-pose interaction
    recall. Missing interaction values on scored poses are a hard error.
    """
    df = master.copy()
    # best (most negative) docking score per ligand across receptors -> ensemble consensus
    df["neg_dock"] = -pd.to_numeric(df["docking_score"], errors="coerce")
    required = {"cnnscore", "key_interaction_recall_consensus"}
    missing = required.difference(df.columns)
    if missing:
        raise ValueError(f"Pose-coupled enrichment requires columns: {sorted(missing)}")
    cnn = pd.to_numeric(df["cnnscore"], errors="coerce")
    recall = pd.to_numeric(df["key_interaction_recall_consensus"], errors="coerce")
    invalid = cnn.notna() & recall.isna()
    if invalid.any():
        raise ValueError(f"{int(invalid.sum())} scored poses lack strict interaction recall")
    df["arm_inter"] = cnn * (1.0 + recall)
    df["arm_gnina"] = cnn
    agg: dict[str, Any] = {
        "neg_dock": "max",
        "label": "first",
        "activity": "first",
        "arm_gnina": "max",
        "arm_inter": "max",
    }
    for col in ["cnnaffinity", "pose_confidence_probability"]:
        if col in df.columns:
            agg[col] = "max"
    per_ligand = df.groupby("lit_pcba_id").agg(agg).reset_index()

    per_ligand["arm_dock"] = per_ligand["neg_dock"]
    if "pose_confidence_probability" in per_ligand.columns:
        per_ligand["arm_pose"] = pd.to_numeric(per_ligand["pose_confidence_probability"], errors="coerce")
    return per_ligand


def step_score(master: pd.DataFrame, out_dir: Path) -> pd.DataFrame:
    per_ligand = build_arms(master)
    per_ligand.to_parquet(out_dir / "enrichment_per_ligand.parquet", index=False)
    arms = [c for c in ["arm_dock", "arm_gnina", "arm_inter", "arm_pose"] if c in per_ligand.columns]
    reports = [enrichment_report(per_ligand, arm) for arm in arms]
    report_df = pd.DataFrame(reports)
    report_df.to_parquet(out_dir / "enrichment_metrics.parquet", index=False)
    report_df.to_csv(out_dir / "enrichment_metrics.csv", index=False)
    (out_dir / "enrichment_metrics.json").write_text(json.dumps(reports, indent=2), encoding="utf-8")
    return report_df


def run_all(work_root: Path, extracted_root: Path, *, max_ligands: int | None = None) -> dict[str, Any]:
    root = project_root()
    tools = default_tools(root)
    work_root = Path(work_root)
    work_root.mkdir(parents=True, exist_ok=True)

    egfr = find_litpcba_egfr(extracted_root)
    if "actives" not in egfr or "decoys" not in egfr:
        raise FileNotFoundError(f"LIT-PCBA EGFR actives/decoys not found under {extracted_root}: {egfr}")

    prepared = step_prepare(egfr["actives"], egfr["decoys"], work_root / "prepared_ligands", tools["obabel"])
    if max_ligands:  # optional cap for a fast validation slice (keep all actives + sample decoys)
        act = prepared[prepared["label"] == 1]
        dec = prepared[prepared["label"] == 0].head(max(0, max_ligands - len(act)))
        prepared = pd.concat([act, dec], ignore_index=True)

    receptors = _receptor_table_with_pdbqt(root)
    master = run_campaign(prepared, receptors, tools, work_root / "campaign")
    report = step_score(master, work_root)
    return {"n_ligands": len(prepared), "n_scored_rows": len(master), "arms": report.to_dict("records")}
