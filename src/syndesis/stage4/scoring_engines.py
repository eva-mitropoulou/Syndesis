from __future__ import annotations

import re
import shutil
import subprocess
from dataclasses import dataclass


@dataclass(frozen=True)
class EngineInfo:
    available: bool
    executable: str
    version: str | None
    warning: str | None = None


def parse_gnina_version(text: str) -> str | None:
    match = re.search(r"gnina\s+v?([0-9][^\s]*)", text, flags=re.IGNORECASE)
    return match.group(0).strip() if match else None


def check_gnina(executable: str, timeout_seconds: int = 60) -> EngineInfo:
    resolved = shutil.which(executable) or executable
    if not resolved:
        return EngineInfo(False, executable, None, "GNINA executable is not configured.")
    try:
        completed = subprocess.run(
            [resolved, "--version"],
            check=False,
            capture_output=True,
            text=True,
            timeout=timeout_seconds,
        )
    except Exception as exc:
        return EngineInfo(False, resolved, None, str(exc))
    text = f"{completed.stdout}\n{completed.stderr}"
    version = parse_gnina_version(text)
    return EngineInfo(completed.returncode == 0, resolved, version, None if completed.returncode == 0 else text[-1000:])

