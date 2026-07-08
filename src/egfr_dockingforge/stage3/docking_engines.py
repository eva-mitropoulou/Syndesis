from __future__ import annotations

import shutil
from dataclasses import dataclass


@dataclass(frozen=True)
class EngineAvailability:
    engine: str
    executable: str | None
    available: bool
    version: str


def check_engine(engine: str, executable: str | None) -> EngineAvailability:
    exe = executable or engine
    found = shutil.which(exe)
    return EngineAvailability(engine=engine, executable=found, available=found is not None, version="available" if found else "unavailable")

