"""Prepare LIT-PCBA (or any SMILES) ligands to docking-ready PDBQT.

SMILES -> RDKit 3D embed (ETKDGv3) + MMFF minimize -> SDF -> OpenBabel PDBQT
(Gasteiger charges), the same PDBQT flavour the Stage-8 screener consumes. Prep is
embarrassingly parallel and CPU-bound, so it runs across a process pool while the
GPU is free; it is idempotent (skips ligands whose PDBQT already exists).

LIT-PCBA ships actives/inactives as SMILES files (``*.smi``: ``SMILES id`` per line);
this module parses those and labels each ligand active/decoy for the enrichment.
"""
from __future__ import annotations

import subprocess
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Any

import pandas as pd


def parse_litpcba_smi(path: Path, label: int) -> list[dict[str, Any]]:
    """Parse a LIT-PCBA .smi file into records. label: 1=active, 0=decoy.
    Lines are 'SMILES id' (whitespace-separated); id may be absent."""
    records = []
    for i, raw in enumerate(path.read_text(encoding="utf-8", errors="replace").splitlines()):
        line = raw.strip()
        if not line:
            continue
        parts = line.split()
        smiles = parts[0]
        mol_id = parts[1] if len(parts) > 1 else f"{path.stem}_{i}"
        records.append({"lit_pcba_id": mol_id, "smiles": smiles, "label": label,
                        "activity": "active" if label == 1 else "decoy"})
    return records


def _prep_one(args: tuple[str, str, str, str]) -> dict[str, Any]:
    """Worker: (ligand_id, smiles, out_pdbqt_path, obabel_path) -> status dict.
    Kept module-level and self-contained so it pickles for ProcessPoolExecutor."""
    ligand_id, smiles, out_pdbqt, obabel = args
    out = Path(out_pdbqt)
    if out.exists() and out.stat().st_size > 0:
        return {"lit_pcba_id": ligand_id, "prepared_pdbqt": str(out), "prep_status": "cached"}
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem
    except Exception as exc:  # noqa: BLE001
        return {"lit_pcba_id": ligand_id, "prepared_pdbqt": None, "prep_status": f"rdkit_import_failed:{exc}"}

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return {"lit_pcba_id": ligand_id, "prepared_pdbqt": None, "prep_status": "invalid_smiles"}
    mol = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = 0xF00D
    if AllChem.EmbedMolecule(mol, params) != 0:
        # retry with random coords for awkward molecules
        params.useRandomCoords = True
        if AllChem.EmbedMolecule(mol, params) != 0:
            return {"lit_pcba_id": ligand_id, "prepared_pdbqt": None, "prep_status": "embed_failed"}
    try:
        AllChem.MMFFOptimizeMolecule(mol, maxIters=1000)
    except Exception:  # noqa: BLE001
        pass  # keep embedded geometry if MMFF unavailable for this molecule
    sdf = out.with_suffix(".sdf")
    out.parent.mkdir(parents=True, exist_ok=True)
    Chem.MolToMolFile(mol, str(sdf))
    try:
        proc = subprocess.run(
            [obabel, str(sdf), "-opdbqt", "-O", str(out), "--partialcharge", "gasteiger"],
            check=False, capture_output=True, text=True, timeout=120,
        )
    except Exception as exc:  # noqa: BLE001
        return {"lit_pcba_id": ligand_id, "prepared_pdbqt": None, "prep_status": f"obabel_error:{exc}"}
    if proc.returncode != 0 or not out.exists() or out.stat().st_size == 0:
        return {"lit_pcba_id": ligand_id, "prepared_pdbqt": None, "prep_status": "obabel_failed"}
    return {"lit_pcba_id": ligand_id, "prepared_pdbqt": str(out), "prep_status": "prepared"}


def prepare_ligands(
    records: list[dict[str, Any]],
    prepared_dir: Path,
    obabel: str,
    *,
    max_workers: int = 10,
) -> pd.DataFrame:
    """Prepare a list of {lit_pcba_id, smiles, label,...} records to PDBQT in
    parallel. Returns the records augmented with prepared_pdbqt + prep_status."""
    prepared_dir.mkdir(parents=True, exist_ok=True)
    jobs = [(r["lit_pcba_id"], r["smiles"], str(prepared_dir / f"{r['lit_pcba_id']}.pdbqt"), obabel) for r in records]
    by_id = {r["lit_pcba_id"]: r for r in records}
    results: list[dict[str, Any]] = []
    with ProcessPoolExecutor(max_workers=max_workers) as pool:
        futs = {pool.submit(_prep_one, j): j[0] for j in jobs}
        for fut in as_completed(futs):
            res = fut.result()
            base = dict(by_id[res["lit_pcba_id"]])
            base.update(res)
            results.append(base)
    return pd.DataFrame(results)
