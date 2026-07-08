from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pandas as pd
import yaml


def project_root() -> Path:
    return Path(__file__).resolve().parents[3]


def resolve_path(path: str | Path, base: Path | None = None) -> Path:
    raw = Path(path).expanduser()
    if raw.is_absolute():
        return raw
    return (base or project_root()) / raw


def ensure_dir(path: str | Path) -> Path:
    resolved = Path(path)
    resolved.mkdir(parents=True, exist_ok=True)
    return resolved


def load_yaml(path: str | Path) -> dict[str, Any]:
    with Path(path).open("r", encoding="utf-8") as handle:
        payload = yaml.safe_load(handle) or {}
    if not isinstance(payload, dict):
        raise TypeError(f"Expected YAML mapping in {path}")
    return payload


def write_json(path: str | Path, payload: dict[str, Any]) -> None:
    target = Path(path)
    ensure_dir(target.parent)
    target.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def write_table(path: str | Path, frame: pd.DataFrame) -> None:
    target = Path(path)
    ensure_dir(target.parent)
    if target.suffix == ".parquet":
        try:
            frame.to_parquet(target, index=False)
        except ImportError as exc:
            raise RuntimeError(
                "Writing parquet requires pyarrow or fastparquet. Install the project "
                "environment from environment.yml or run `pip install pyarrow` in your "
                "active environment."
            ) from exc
    elif target.suffix == ".csv":
        frame.to_csv(target, index=False)
    else:
        raise ValueError(f"Unsupported table extension: {target.suffix}")
