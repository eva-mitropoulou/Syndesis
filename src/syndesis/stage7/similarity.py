from __future__ import annotations

import numpy as np
from rdkit import Chem, DataStructs
from rdkit.Chem import rdFingerprintGenerator


MORGAN = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)


def fingerprint(smiles: str):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    return MORGAN.GetFingerprint(mol)


def prepare_known(known: list[tuple[str, str, str, float | None, str | None]]) -> list[tuple[str, str, object, float | None, str | None]]:
    prepared = []
    for molecule_id, source, known_smiles, activity, endpoint in known:
        fp = fingerprint(known_smiles)
        if fp is not None:
            prepared.append((molecule_id, source, fp, activity, endpoint))
    return prepared


def closest(smiles: str, known: list[tuple[str, str, object, float | None, str | None]]) -> tuple[str, str, float | None, str | None, float]:
    fp = fingerprint(smiles)
    if fp is None or not known:
        return "", "", None, None, np.nan
    best = ("", "", None, None, -1.0)
    for molecule_id, source, kfp, activity, endpoint in known:
        sim = float(DataStructs.TanimotoSimilarity(fp, kfp))
        if sim > best[4]:
            best = (molecule_id, source, activity, endpoint, sim)
    return best


def novelty_bucket(similarity: float, same_key: bool, scaffold_match: bool) -> str:
    if same_key:
        return "known_duplicate"
    if similarity >= 0.70:
        return "close_analog"
    if 0.40 <= similarity < 0.70:
        return "medium_similarity"
    if scaffold_match:
        return "low_similarity_known_scaffold"
    return "scaffold_novel"
