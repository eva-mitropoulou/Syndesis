"""Audit CDK2 native-ligand chemical overlap with DUD-E actives."""
from __future__ import annotations

import urllib.request
from pathlib import Path

import pandas as pd
from rdkit import Chem
from rdkit.Chem import rdFingerprintGenerator
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit.DataStructs import BulkTanimotoSimilarity

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results_showcase/submission_robustness/cdk2_native_dude_overlap_audit.csv"
BASE = ROOT / "data/external/cdk2"
CCD_ROOT = ROOT / "data/references/cdk2_native_ligands"


def ccd_molecule(component_id: str) -> tuple[Chem.Mol, Path, str]:
    CCD_ROOT.mkdir(parents=True, exist_ok=True)
    url = f"https://files.rcsb.org/ligands/download/{component_id}_ideal.sdf"
    path = CCD_ROOT / f"{component_id}_ideal.sdf"
    if not path.exists():
        with urllib.request.urlopen(url, timeout=30) as response:  # noqa: S310
            path.write_bytes(response.read())
    molecule = Chem.SDMolSupplier(str(path), removeHs=True)[0]
    if molecule is None:
        raise RuntimeError(f"RCSB CCD ideal structure is not a valid molecule: {path}")
    return molecule, path, url


def main() -> int:
    active_rows = []
    with (BASE / "actives_final.ism").open() as handle:
        for line in handle:
            fields = line.split()
            if len(fields) >= 2:
                molecule = Chem.MolFromSmiles(fields[0])
                if molecule is not None:
                    active_rows.append((fields[1], molecule))
    generator = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    active_fingerprints = [generator.GetFingerprint(molecule) for _, molecule in active_rows]
    active_scaffolds = [MurckoScaffold.MurckoScaffoldSmiles(mol=molecule) for _, molecule in active_rows]

    ensemble = pd.read_parquet(BASE / "cdk2_receptor_ensemble.parquet")
    ensemble = ensemble[ensemble["selected_flag"].astype(bool)]
    rows = []
    for record in ensemble.to_dict("records"):
        receptor_id = record["receptor_id"]
        molecule, path, url = ccd_molecule(record["ligand_comp_id"])
        smiles = Chem.MolToSmiles(molecule, isomericSmiles=True)
        fingerprint = generator.GetFingerprint(molecule)
        similarities = BulkTanimotoSimilarity(fingerprint, active_fingerprints)
        nearest_index = max(range(len(similarities)), key=similarities.__getitem__)
        scaffold = MurckoScaffold.MurckoScaffoldSmiles(mol=molecule)
        rows.append({
            "target": "CDK2",
            "receptor_id": receptor_id,
            "pdb_id": record["pdb_id"],
            "native_ligand": record["ligand_comp_id"],
            "native_smiles": smiles,
            "native_structure_path": str(path.relative_to(ROOT)),
            "native_structure_source_url": url,
            "n_dude_actives": len(active_rows),
            "nearest_active_id": active_rows[nearest_index][0],
            "max_ecfp4_tanimoto": similarities[nearest_index],
            "mean_ecfp4_tanimoto": sum(similarities) / len(similarities),
            "actives_ge_0_7": sum(value >= 0.7 for value in similarities),
            "actives_ge_0_5": sum(value >= 0.5 for value in similarities),
            "native_bemis_murcko_scaffold": scaffold,
            "scaffold_matches_an_active": scaffold in active_scaffolds,
            "fingerprint": "Morgan/ECFP4 radius=2, 2048 bits",
        })
    OUT.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(OUT, index=False)
    print(pd.DataFrame(rows).to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
