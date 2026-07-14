from __future__ import annotations

import pandas as pd
from rdkit import Chem
from rdkit.Chem import BRICS


def brics_fragments(smiles: str) -> list[str]:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return []
    return sorted(BRICS.BRICSDecompose(mol))
