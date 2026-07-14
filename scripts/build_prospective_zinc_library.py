"""Rebuild and verify the deterministic ZINC prospective screening library."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import pandas as pd
from rdkit import Chem, RDLogger, rdBase
from rdkit.Chem import Crippen, Descriptors, Lipinski, QED, rdMolDescriptors
from rdkit.Chem.Scaffolds import MurckoScaffold

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_SOURCE = ROOT / "data/external/prospective/zinc_tranches"
ALLOWED_ELEMENTS = {1, 5, 6, 7, 8, 9, 14, 15, 16, 17, 35, 53}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def descriptors(smiles: str) -> dict | None:
    molecule = Chem.MolFromSmiles(smiles)
    if molecule is None:
        return None
    try:
        return {
            "mw": Descriptors.MolWt(molecule),
            "logp": Crippen.MolLogP(molecule),
            "hbd": Lipinski.NumHDonors(molecule),
            "hba": Lipinski.NumHAcceptors(molecule),
            "rot": Lipinski.NumRotatableBonds(molecule),
            "tpsa": rdMolDescriptors.CalcTPSA(molecule),
            "arom_n": sum(atom.GetIsAromatic() and atom.GetAtomicNum() == 7 for atom in molecule.GetAtoms()),
            "qed": QED.qed(molecule),
            "scaffold": MurckoScaffold.MurckoScaffoldSmiles(mol=molecule),
            "single_fragment": "." not in smiles,
            "allowed_elements": all(atom.GetAtomicNum() in ALLOWED_ELEMENTS for atom in molecule.GetAtoms()),
        }
    except RuntimeError:
        return None


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source-root", type=Path, default=DEFAULT_SOURCE)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--verify-existing", type=Path)
    args = parser.parse_args()
    RDLogger.DisableLog("rdApp.error")

    source_files = sorted(args.source_root.glob("tr_*.smi"))
    if len(source_files) != 36:
        raise RuntimeError(f"Expected 36 ZINC tranche files, found {len(source_files)}")
    manifest_rows = []
    frames = []
    for path in source_files:
        frame = pd.read_csv(path, sep=r"\s+", dtype=str)
        frame["source_file"] = path.name
        frames.append(frame)
        manifest_rows.append({
            "source_file": path.name,
            "tranche": path.stem.removeprefix("tr_"),
            "n_records": len(frame),
            "bytes": path.stat().st_size,
            "sha256": sha256(path),
            "retrieval_date": "2026-07-14",
            "source": "ZINC22/CartBlanche22 tranche export",
            "source_url": "https://cartblanche22.docking.org/",
        })
    raw = pd.concat(frames, ignore_index=True)
    n_raw = len(raw)
    raw = raw.drop_duplicates("zinc_id", keep="first").copy()
    values = raw["smiles"].map(descriptors)
    valid = values.notna()
    calculated = pd.DataFrame(values[valid].tolist(), index=raw.index[valid])
    raw = raw.join(calculated)
    filtered = raw[
        valid
        & raw["single_fragment"].eq(True)
        & raw["allowed_elements"].eq(True)
        & raw["mw"].between(300.0, 500.0)
        & raw["logp"].between(1.0, 5.0)
        & raw["hbd"].le(3)
        & raw["hba"].le(8)
        & raw["rot"].le(8)
        & raw["tpsa"].le(110.0)
        & raw["arom_n"].ge(1)
        & raw["qed"].ge(0.4)
    ].copy()
    n_property_pass = len(filtered)
    filtered = filtered[filtered.groupby("scaffold", sort=False).cumcount() < 2].copy()
    n_diversity_pass = len(filtered)
    library = filtered.head(2000).copy()
    library["mw"] = library["mw"].round(1)
    library["logp"] = library["logp"].round(2)
    library["tpsa"] = library["tpsa"].round(1)
    library["qed"] = library["qed"].round(3)
    library["lit_pcba_id"] = library["zinc_id"]
    library["label"] = -1
    library["activity"] = "prospective"
    columns = [
        "smiles", "zinc_id", "mw", "logp", "hbd", "hba", "rot", "tpsa",
        "arom_n", "qed", "scaffold", "lit_pcba_id", "label", "activity",
    ]
    library = library[columns].reset_index(drop=True)

    manifest = pd.DataFrame(manifest_rows)
    manifest_path = ROOT / "data/references/prospective_zinc_source_manifest.csv"
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    manifest.to_csv(manifest_path, index=False)
    audit = {
        "n_source_files": len(source_files),
        "n_raw_records": n_raw,
        "n_unique_zinc_ids": len(raw),
        "n_property_filter_pass": n_property_pass,
        "n_scaffold_cap_pass": n_diversity_pass,
        "max_molecules_per_scaffold": 2,
        "n_selected": len(library),
        "selection_order": "source filename lexical order, then source row order",
        "random_sampling": False,
        "chembl_exclusion_applied": False,
        "rdkit_version": rdBase.rdkitVersion,
    }

    if args.verify_existing:
        existing = pd.read_parquet(args.verify_existing) if args.verify_existing.suffix == ".parquet" else pd.read_csv(args.verify_existing)
        identity_columns = ["smiles", "zinc_id", "scaffold", "lit_pcba_id", "label", "activity"]
        pd.testing.assert_frame_equal(
            library[identity_columns], existing[identity_columns], check_dtype=False
        )
        tolerances = {"mw": 0.11, "logp": 0.011, "tpsa": 0.11, "qed": 0.0011}
        drift = {}
        for column, tolerance in tolerances.items():
            difference = (library[column] - existing[column]).abs()
            if difference.max() > tolerance:
                raise AssertionError(
                    f"{column} descriptor drift {difference.max()} exceeds tolerance {tolerance}"
                )
            drift[column] = {
                "n_different": int(difference.gt(0).sum()),
                "maximum_absolute_difference": float(difference.max()),
            }
        for column in ["hbd", "hba", "rot", "arom_n"]:
            pd.testing.assert_series_equal(library[column], existing[column], check_dtype=False)
        audit["verification"] = {
            "existing_library": str(args.verify_existing),
            "exact_identity_and_order_match": True,
            "descriptor_drift": drift,
        }
    audit_path = ROOT / "results_showcase/submission_robustness/prospective_library_rebuild_audit.json"
    audit_path.parent.mkdir(parents=True, exist_ok=True)
    audit_path.write_text(json.dumps(audit, indent=2) + "\n")
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        if args.output.suffix == ".parquet":
            library.to_parquet(args.output, index=False)
        else:
            library.to_csv(args.output, index=False)
    print(json.dumps(audit), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
