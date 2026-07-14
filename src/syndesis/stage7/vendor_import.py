from __future__ import annotations

from io import StringIO
from pathlib import Path
from typing import Any

import pandas as pd
import requests


def fetch_cartblanche_similarity(smiles: str, distance: int, timeout: int = 120) -> pd.DataFrame:
    response = requests.get(
        "https://cartblanche22.docking.org/smiles.txt",
        files={
            "smiles": (None, smiles),
            "dist": (None, str(distance)),
            "output_fields": (None, "zinc_id,smiles,tranche,catalogs"),
        },
        timeout=timeout,
    )
    response.raise_for_status()
    text = response.text.strip()
    if text.startswith("{") and "task" in text:
        raise RuntimeError(f"CartBlanche returned asynchronous task instead of table: {text}")
    if not text or text.startswith("<!doctype"):
        raise RuntimeError("CartBlanche did not return a parseable vendor table")
    return pd.read_csv(StringIO(text), sep="\t")


def read_vendor_file(path: str | Path) -> pd.DataFrame:
    p = Path(path)
    if p.suffix == ".parquet":
        return pd.read_parquet(p)
    if p.suffix in {".csv", ".tsv"}:
        return pd.read_csv(p, sep="\t" if p.suffix == ".tsv" else ",")
    return pd.read_csv(p, sep=None, engine="python", names=["smiles", "vendor_catalog_id"])
