"""Package path-independent MD topology and parameter files for the paper archive."""
from __future__ import annotations

import argparse
import hashlib
import re
import shutil
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
SYSTEM_FILES = [
    "complex.gro",
    "posre.itp",
    "ions.mdp",
    "minimization.mdp",
    "nvt_equilibration.mdp",
    "npt_equilibration.mdp",
    "production_quick.mdp",
    "production_quick_rep02.mdp",
    "production_quick_rep03.mdp",
]
PARAMETER_SUFFIXES = {
    "_GMX.itp": "ligand_GMX.itp",
    "_GMX.gro": "ligand_GMX.gro",
    ".mol2": "ligand_gaff2_am1bcc.mol2",
    "_AC.frcmod": "ligand_AC.frcmod",
    "_AC.prmtop": "ligand_AC.prmtop",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def sanitize_topology(source: Path, target: Path) -> None:
    text = source.read_text()
    marker = "; Include forcefield parameters"
    if marker not in text:
        raise RuntimeError(f"Topology lacks expected standalone marker: {source}")
    text = text[text.index(marker):]
    text, replacements = re.subn(r'#include "[^"]+_GMX\.itp"', '#include "ligand_GMX.itp"', text, count=1)
    if replacements != 1:
        raise RuntimeError(f"Topology does not contain exactly one GAFF ligand include: {source}")
    text, replacements = re.subn(r'#include "[^"]+/posre\.itp"', '#include "posre.itp"', text, count=1)
    if replacements != 1:
        raise RuntimeError(f"Topology does not contain exactly one protein-restraint include: {source}")
    target.write_text(text)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("source_root", type=Path, help="Original egfr_md_work directory")
    parser.add_argument(
        "--output-root",
        type=Path,
        default=ROOT / "results" / "md" / "reproducibility",
    )
    args = parser.parse_args()
    manifest = pd.read_csv(ROOT / "results" / "md" / "md_candidate_manifest.csv")
    args.output_root.mkdir(parents=True, exist_ok=True)
    checksum_rows = []
    for row in manifest.itertuples(index=False):
        system_id = f"mdsys_{row.md_candidate_id}"
        system_source = args.source_root / system_id
        parameter_source = args.source_root / "ligand_parameters" / row.molecule_id / f"{row.molecule_id}.acpype"
        system_target = args.output_root / row.md_candidate_id
        system_target.mkdir(parents=True, exist_ok=True)
        sanitize_topology(system_source / "topol.top", system_target / "topol.top")
        for filename in SYSTEM_FILES:
            shutil.copy2(system_source / filename, system_target / filename)
        for suffix, target_name in PARAMETER_SUFFIXES.items():
            shutil.copy2(parameter_source / f"{row.molecule_id}{suffix}", system_target / target_name)
        log_text = (parameter_source / "acpype.log").read_text(errors="replace")
        log_text = re.sub(
            r"/home/[^\s]+/[A-Za-z0-9_.-]+",
            "${PROJECT_ROOT}",
            log_text,
        )
        log_text = log_text.replace(str(args.source_root), "${MD_WORK_ROOT}")
        (system_target / "acpype.log").write_text(log_text)
        shutil.copy2(parameter_source / "sqm.out", system_target / "sqm.out")
        for path in sorted(system_target.iterdir()):
            checksum_rows.append({
                "md_candidate_id": row.md_candidate_id,
                "molecule_id": row.molecule_id,
                "relative_path": str(path.relative_to(args.output_root)),
                "sha256": sha256(path),
                "size_bytes": path.stat().st_size,
            })
    pd.DataFrame(checksum_rows).to_csv(args.output_root / "checksums.csv", index=False)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
