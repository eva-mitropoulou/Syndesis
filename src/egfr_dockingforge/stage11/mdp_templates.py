from __future__ import annotations

from pathlib import Path
from typing import Any


def minimization_mdp(config: dict[str, Any]) -> str:
    return f"""integrator = steep
emtol = 1000.0
emstep = 0.01
nsteps = {int(config['md']['minimization_max_steps'])}
coulombtype = PME
rcoulomb = 1.2
rvdw = 1.2
vdwtype = switch
rvdw-switch = 1.0
constraints = none
tcoupl = no
pcoupl = no
"""


def nvt_mdp(config: dict[str, Any]) -> str:
    steps = int(float(config["md"]["nvt_ps"]) / (float(config["md"]["timestep_fs"]) / 1000.0))
    tc_grps = config["md"].get("thermostat_groups", "Protein Non-Protein")
    return f"""integrator = md
define = -DPOSRES
nsteps = {steps}
dt = {float(config['md']['timestep_fs'])/1000.0}
continuation = no
constraint_algorithm = lincs
constraints = h-bonds
tcoupl = V-rescale
; NOTE: for a rigorous protocol the ligand should be grouped with the protein via a custom index group (e.g. Protein_UNL); the default lumps the ligand with solvent/ions.
tc-grps = {tc_grps}
tau_t = 0.1 0.1
ref_t = {config['md']['temperature_k']} {config['md']['temperature_k']}
pcoupl = no
coulombtype = PME
rcoulomb = 1.2
rvdw = 1.2
gen_vel = yes
"""


def npt_mdp(config: dict[str, Any]) -> str:
    steps = int(float(config["md"]["npt_ps"]) / (float(config["md"]["timestep_fs"]) / 1000.0))
    tc_grps = config["md"].get("thermostat_groups", "Protein Non-Protein")
    return f"""integrator = md
define = -DPOSRES
nsteps = {steps}
dt = {float(config['md']['timestep_fs'])/1000.0}
continuation = yes
constraint_algorithm = lincs
constraints = h-bonds
tcoupl = V-rescale
; NOTE: for a rigorous protocol the ligand should be grouped with the protein via a custom index group (e.g. Protein_UNL); the default lumps the ligand with solvent/ions.
tc-grps = {tc_grps}
tau_t = 0.1 0.1
ref_t = {config['md']['temperature_k']} {config['md']['temperature_k']}
pcoupl = Parrinello-Rahman
pcoupltype = isotropic
tau_p = 2.0
ref_p = {config['md']['pressure_bar']}
compressibility = 4.5e-5
coulombtype = PME
rcoulomb = 1.2
rvdw = 1.2
gen_vel = no
"""


def production_mdp(config: dict[str, Any], replicate: bool = False) -> str:
    ns = float(config["md"]["replicate_production_ns"] if replicate else config["md"]["quick_production_ns"])
    steps = int(ns * 1000 / (float(config["md"]["timestep_fs"]) / 1000.0))
    tc_grps = config["md"].get("thermostat_groups", "Protein Non-Protein")
    return f"""integrator = md
nsteps = {steps}
dt = {float(config['md']['timestep_fs'])/1000.0}
continuation = yes
constraint_algorithm = lincs
constraints = h-bonds
tcoupl = V-rescale
; NOTE: for a rigorous protocol the ligand should be grouped with the protein via a custom index group (e.g. Protein_UNL); the default lumps the ligand with solvent/ions.
tc-grps = {tc_grps}
tau_t = 0.1 0.1
ref_t = {config['md']['temperature_k']} {config['md']['temperature_k']}
pcoupl = Parrinello-Rahman
pcoupltype = isotropic
tau_p = 2.0
ref_p = {config['md']['pressure_bar']}
compressibility = 4.5e-5
coulombtype = PME
rcoulomb = 1.2
rvdw = 1.2
compressed-x-grps = System
nstxout-compressed = 5000
gen_vel = no
; no ligand restraints in production
"""


def ions_mdp(config: dict[str, Any]) -> str:
    return """integrator = steep
emtol = 1000.0
emstep = 0.01
nsteps = 500
coulombtype = PME
rcoulomb = 1.2
rvdw = 1.2
constraints = none
tcoupl = no
pcoupl = no
"""


def write_mdp_templates(config: dict[str, Any], md_dir: Path) -> dict[str, Path]:
    md_dir.mkdir(parents=True, exist_ok=True)
    files = {
        "minimization": (md_dir / "minimization.mdp", minimization_mdp(config)),
        "nvt_equilibration": (md_dir / "nvt_equilibration.mdp", nvt_mdp(config)),
        "npt_equilibration": (md_dir / "npt_equilibration.mdp", npt_mdp(config)),
        "production_quick": (md_dir / "production_quick.mdp", production_mdp(config, False)),
        "production_replicate": (md_dir / "production_replicate.mdp", production_mdp(config, True)),
        "ions": (md_dir / "ions.mdp", ions_mdp(config)),
    }
    return {name: _write(path, text) for name, (path, text) in files.items()}


def _write(path: Path, text: str) -> Path:
    path.write_text(text, encoding="utf-8")
    return path
