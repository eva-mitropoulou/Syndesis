"""Prep a CDK2 receptor ensemble for the second-target enrichment spine.

Cross-family generalization test (CDK2 = Ser/Thr kinase vs EGFR tyrosine kinase).
Builds the same receptor-table contract the enrichment orchestrator consumes:
  receptor_id, selected_flag, suggested_docking_box_center/size, docking_format_file.

For each PDB in the ensemble we:
  1. split chain A ATOM records (drop cyclin chains, waters, other hetero) -> protein PDB.
     This represents 1QMZ without HETATM TPO160; it is not a phosphothreonine model.
  2. extract the named ATP-site ligand (by 3-letter resname) -> ligand PDB + SDF
  3. obabel protein PDB -> rigid receptor PDBQT (same recipe as stage3 docking_prep:
     -xr -p 7.4 --partialcharge gasteiger)
  4. box center = ligand heavy-atom centroid; box size = per-axis (ligand extent + 2*pad),
     floored at 18 A (matches EGFR ensemble box sizing)

Ensemble (holo, ATP-site ligand, chain A), spanning CDK2's activation axis:
  1QMZ active/aC-in/cyclinA/pT160 (ATP) | 1FIN cyclin-bound unphospho (ATP)
  2A4L monomeric roscovitine (RRC)      | 1AQ1 monomeric staurosporine (STU)
  1PXN monomeric ATP-site inhibitor (CK6)
Held out as DUD-E reference (NOT in ensemble): 1H00.

Pure CPU. Writes under /mnt/e/cdk2_dude/receptors + a cdk2_receptor_ensemble.parquet.
"""
from __future__ import annotations

import ast
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

PROJ = Path("/home/dimit/eva/Computational-Chemistry/protein_docking/syndesis")
sys.path.insert(0, str(PROJ / "src"))
from syndesis.enrichment.run_enrichment import default_tools  # noqa: E402

BASE = Path("/mnt/e/cdk2_dude")
ENS_PDB = BASE / "ensemble"
RECEPTORS = BASE / "receptors"
OBABEL = default_tools(PROJ)["obabel"]  # canonical obabel (egfr-cadd env)
PAD = 8.0        # Angstrom padding each side of the ligand extent
MIN_BOX = 18.0   # floor per axis (matches EGFR ensemble)

# (pdb_id, chain, ligand_resname, state_label, role, selected)
ENSEMBLE = [
    ("1QMZ", "A", "ATP", "active_like",   "active_cyclin_pT160", True),
    ("1FIN", "A", "ATP", "active_like",   "cyclin_unphospho",    True),
    ("2A4L", "A", "RRC", "inactive_like", "monomeric_roscovitine", True),
    ("1AQ1", "A", "STU", "inactive_like", "monomeric_staurosporine", True),
    ("1PXN", "A", "CK6", "inactive_like", "monomeric_atp_inhibitor", True),
    ("1H00", "A", "FCP", "reference",     "dude_reference_holdout", False),
]


def _read_pdb_lines(pdb: Path) -> list[str]:
    return pdb.read_text(errors="replace").splitlines()


def _select_altlocs(records: list[str]) -> list[str]:
    """Keep one deposited alternate location per atom, preferring occupancy.

    Open Babel must not receive both conformers for a protein atom: it can retain
    duplicate coordinates while losing the residue assignment in PDBQT. The chosen
    conformer is written with a blank alternate-location field so downstream tools
    see one chemically coherent receptor model.
    """
    selected: dict[tuple[str, str, str, str, str], tuple[float, int, str]] = {}
    for index, line in enumerate(records):
        key = (line[21], line[22:26], line[26], line[17:20], line[12:16])
        altloc = line[16]
        try:
            occupancy = float(line[54:60])
        except ValueError:
            occupancy = 0.0
        # Prefer blank, then A, only when occupancies are equal.
        altloc_priority = 2 if altloc == " " else 1 if altloc == "A" else 0
        candidate = (occupancy, altloc_priority, index)
        current = selected.get(key)
        if current is None or candidate > current[:3]:
            selected[key] = (occupancy, altloc_priority, index)
    return [records[selected[key][2]][:16] + " " + records[selected[key][2]][17:] for key in sorted(selected, key=lambda key: selected[key][2])]


def split_protein_and_ligand(pdb_id: str, chain: str, lig_resn: str) -> dict:
    low = pdb_id.lower()
    src = ENS_PDB / f"{low}.pdb"
    lines = _read_pdb_lines(src)
    prot_lines, lig_lines = [], []
    lig_coords = []
    for ln in lines:
        rec = ln[:6].strip()
        if rec == "ATOM":
            # protein backbone/side-chain atoms of the requested chain only
            if ln[21] == chain:
                prot_lines.append(ln)
        elif rec == "HETATM":
            resn = ln[17:20].strip()
            ch = ln[21]
            if resn == lig_resn and ch == chain:
                lig_lines.append(ln)
                try:
                    lig_coords.append((float(ln[30:38]), float(ln[38:46]), float(ln[46:54])))
                except ValueError:
                    pass
    prot_lines = _select_altlocs(prot_lines)
    if lig_lines:
        lig_lines = _select_altlocs(lig_lines)
        lig_coords = [
            (float(line[30:38]), float(line[38:46]), float(line[46:54]))
            for line in lig_lines
        ]
    if not lig_coords:
        # some ligands are modelled on a different chain id; retry ignoring chain for HETATM
        for ln in lines:
            if ln[:6].strip() == "HETATM" and ln[17:20].strip() == lig_resn:
                lig_lines.append(ln)
                lig_coords.append((float(ln[30:38]), float(ln[38:46]), float(ln[46:54])))
    if not prot_lines:
        raise RuntimeError(f"{pdb_id}: no ATOM lines for chain {chain}")
    if not lig_coords:
        raise RuntimeError(f"{pdb_id}: ligand {lig_resn} not found")

    RECEPTORS.mkdir(parents=True, exist_ok=True)
    prot_pdb = RECEPTORS / f"{low}_{chain}_protein.pdb"
    prot_pdb.write_text("\n".join(prot_lines) + "\nEND\n")
    lig_pdb = RECEPTORS / f"{low}_{chain}_{lig_resn}_ligand.pdb"
    lig_pdb.write_text("\n".join(lig_lines) + "\nEND\n")

    coords = np.array(lig_coords, dtype=float)
    center = coords.mean(axis=0)
    extent = coords.max(axis=0) - coords.min(axis=0)
    size = np.maximum(extent + 2 * PAD, MIN_BOX)
    return {
        "protein_pdb": prot_pdb, "ligand_pdb": lig_pdb,
        "center": [round(float(x), 3) for x in center],
        "size": [round(float(x), 3) for x in size],
        "n_lig_atoms": len(lig_coords),
    }


def make_receptor_pdbqt(protein_pdb: Path, receptor_id: str) -> Path:
    out = RECEPTORS / f"{receptor_id}.pdbqt"
    # same recipe as stage3/docking_prep.prepare_receptors
    cmd = [OBABEL, str(protein_pdb), "-O", str(out), "-xr", "-p", "7.4",
           "--partialcharge", "gasteiger"]
    res = subprocess.run(cmd, capture_output=True, text=True)
    if not out.exists() or out.stat().st_size == 0:
        raise RuntimeError(f"obabel receptor prep failed for {receptor_id}: {res.stderr[-400:]}")
    return out


def make_ligand_sdf(ligand_pdb: Path, receptor_id: str) -> Path | None:
    """Native ligand -> SDF (coords preserved, protonated) for the interaction consensus."""
    out = RECEPTORS / f"{receptor_id}_native_ligand.sdf"
    cmd = [OBABEL, str(ligand_pdb), "-O", str(out), "-p", "7.4"]
    subprocess.run(cmd, capture_output=True, text=True)
    return out if out.exists() and out.stat().st_size > 0 else None


def main() -> int:
    rows = []
    for pdb_id, chain, lig, state, role, selected in ENSEMBLE:
        rid = f"{pdb_id.lower()}_{chain.lower()}_{lig.lower()}"
        try:
            sp = split_protein_and_ligand(pdb_id, chain, lig)
            pdbqt = make_receptor_pdbqt(sp["protein_pdb"], rid)
            sdf = make_ligand_sdf(sp["ligand_pdb"], rid)
            n_rec_atoms = sum(1 for _ in pdbqt.read_text(errors="replace").splitlines()
                              if _.startswith(("ATOM", "HETATM")))
            rows.append({
                "receptor_id": rid, "pdb_id": pdb_id, "auth_asym_id": chain,
                "ligand_comp_id": lig, "state_stratum": state, "selected_role": role,
                "modified_residue_handling": "TPO160_omitted_by_atom_record_filter" if pdb_id == "1QMZ" else "standard_atom_records_only",
                "selected_flag": bool(selected),
                "suggested_docking_box_center": sp["center"],
                "suggested_docking_box_size": sp["size"],
                "docking_format_file": str(pdbqt),
                "native_ligand_sdf_path": str(sdf) if sdf else None,
                "n_ligand_atoms": sp["n_lig_atoms"], "n_receptor_atoms": n_rec_atoms,
            })
            print(f"  OK {rid}: box c={sp['center']} s={sp['size']} "
                  f"rec_atoms={n_rec_atoms} lig_atoms={sp['n_lig_atoms']} "
                  f"selected={selected}", flush=True)
        except Exception as exc:  # noqa: BLE001
            print(f"  FAIL {pdb_id} ({lig}): {exc}", flush=True)
    if not rows:
        print("\nNO receptors prepared — aborting (check obabel path / ligand resnames)")
        return 1
    df = pd.DataFrame(rows)
    out = BASE / "cdk2_receptor_ensemble.parquet"
    df.to_parquet(out, index=False)
    df.to_csv(BASE / "cdk2_receptor_ensemble.csv", index=False)
    n_sel = int(df["selected_flag"].sum()) if "selected_flag" in df else 0
    print(f"\nwrote {out}  ({n_sel} selected / {len(df)} total)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
