from __future__ import annotations

import shutil


def converter_available(executable: str = "cgenff_charmm2gmx") -> bool:
    return shutil.which(executable) is not None
