"""Rerun the CDK2 1PXN arm after deterministic alternate-location selection.

The deposited 1PXN Lys33 side chain has A/B alternate locations. The original
PDBQT retained both, whereas the ProLIF receptor retained only the 0.98-occupancy
A conformer. This driver reruns only 1PXN using the corrected single-conformer
receptor and writes a separate, fail-closed campaign checkpoint.
"""
from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path

import pandas as pd

from syndesis.enrichment.orchestrator import _box_for_receptor, run_receptor
from syndesis.enrichment.run_enrichment import default_tools


ROOT = Path(__file__).resolve().parents[1]
CDK2 = Path(os.environ.get("SYNTHESIS_CDK2_WORK", "cdk2_dude"))
OLD_WORK = Path(os.environ.get("SYNTHESIS_CDK2_BASELINE", "cdk2_enrichment"))
WORK = Path(os.environ.get("SYNTHESIS_1PXN_RERUN", "cdk2_1pxn_altloc_rerun"))
RECEPTOR_ID = "1pxn_a_ck6"


def _heavy_atom_count(path: Path) -> int:
    count = 0
    for line in path.read_text(errors="replace").splitlines():
        if not line.startswith(("ATOM", "HETATM")):
            continue
        name = line[12:16].strip()
        element = line[76:78].strip().upper()
        if element == "H" or (not element and name.startswith("H")):
            continue
        count += 1
    return count


def main() -> int:
    tools = default_tools(ROOT)
    receptor = CDK2 / "receptors" / f"{RECEPTOR_ID}.pdbqt"
    prolif = CDK2 / "stage5" / "prolif_proteins" / "receptors__1pxn_a_protein.h.pdb"
    if not receptor.exists() or not prolif.exists():
        raise FileNotFoundError("Corrected 1PXN receptor or its ProLIF model is missing")
    if _heavy_atom_count(receptor) != _heavy_atom_count(prolif):
        raise RuntimeError("1PXN docking and ProLIF heavy-atom counts differ; rerun blocked")

    prepared = pd.read_parquet(OLD_WORK / "ligand_prep.parquet")
    # The source Lit-PCBA table contains 19 duplicate molecule IDs.  The frozen
    # CDK2 benchmark retains one prepared structure per distinct molecule.
    ligands = (
        prepared[prepared["prep_status"].isin(["prepared", "cached"])]
        .drop_duplicates("lit_pcba_id", keep="first")
        .copy()
    )
    if len(ligands) != 28_296:
        raise RuntimeError(f"Expected 28,296 prepared CDK2 ligands, found {len(ligands)}")
    receptors = pd.read_parquet(CDK2 / "cdk2_receptor_ensemble.parquet")
    row = receptors[receptors["receptor_id"].eq(RECEPTOR_ID)]
    if len(row) != 1:
        raise RuntimeError("1PXN receptor row missing or duplicated")

    WORK.mkdir(parents=True, exist_ok=True)
    manifest = {
        "receptor_id": RECEPTOR_ID,
        "reason": "remove 0.02-occupancy Lys33-B alternate location from docking receptor",
        "docking_receptor": str(receptor),
        "docking_receptor_sha256": hashlib.sha256(receptor.read_bytes()).hexdigest(),
        "docking_heavy_atoms": _heavy_atom_count(receptor),
        "prolif_heavy_atoms": _heavy_atom_count(prolif),
        "n_ligands": len(ligands),
        "unidock_num_modes": 9,
        "unidock_seed": 807,
        "gnina_image": tools["gnina_image"],
    }
    (WORK / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")

    table = run_receptor(
        RECEPTOR_ID,
        str(receptor),
        _box_for_receptor(receptors, RECEPTOR_ID),
        ligands,
        tools,
        WORK / "campaign",
        num_modes=9,
        seed=807,
    )
    table.to_parquet(WORK / f"scores_{RECEPTOR_ID}.parquet", index=False)
    table.to_csv(WORK / f"scores_{RECEPTOR_ID}.csv", index=False)
    if table["cnnscore"].notna().sum() == 0:
        raise RuntimeError("1PXN rerun produced no CNNscores")
    (WORK / "RERUN_DONE").write_text("done\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
