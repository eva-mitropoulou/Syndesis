# Reproducibility Notes

Practical setup notes for reproducing the Syndesis pipeline. See
`environment.yml`, `pyproject.toml`, and `configs/tools.example.yaml` for the
canonical dependency lists.

## (a) Environment

The pipeline requires a full scientific Python stack, not just the base CLI
deps. Preferred install is conda/micromamba from `environment.yml`, then an
editable install of this package:

```
micromamba env create -f environment.yml
micromamba activate syndesis
pip install -e .
```

Key Python packages that MUST be importable: `rdkit`, `MDAnalysis`,
`scikit-learn` (`sklearn`), `lightgbm`, `xgboost`, `catboost`, `prolif`,
`scipy`, `shap`, `openpyxl`, plus `openmm` and `pdbfixer` (Stage 5 ProLIF
engine). `openmm`/`pdbfixer` are conda-forge packages and are not reliably
pip-installable; install them via conda. `catboost` is a real dependency
(Stage 6 ranker) even if not present in every runtime; it is declared in
`pyproject.toml`.

## (b) GROMACS 2026 TPR + MDAnalysis

MDAnalysis cannot parse GROMACS 2026 TPR files (tpx version 138). Stage 11
trajectory analysis therefore takes its topology from the `.gro` file and
delegates PBC correction and fitting to `gmx trjconv` (`-pbc mol -center`
then `-fit rot+trans`) before loading the corrected `.xtc` into MDAnalysis.
The TPR is still needed as the `-s` reference for `trjconv`.

## (c) Deployment-specific tool paths in configs

The following config keys hold absolute paths from the original author's
machine and MUST be edited per deployment (marked `# EDIT ME:` in-file). Copy
`configs/tools.example.yaml` to `configs/tools.yaml` as a template.

- `configs/stage3_redocking_crossdocking.yaml`: `docking.engine_executable_path.unidock`, `prep.obabel_path`
- `configs/stage4_rescoring.yaml`: `gnina.executable`
- `configs/stage8_candidate_screening.yaml`: `docking.executable` (unidock), `docking.obabel`, `gnina.executable`
- `configs/stage11_md_stress_test.yaml`: `forcefield.gromacs_executable` and the AmberTools/CGenFF tool paths

The optional remote CPU worker in `configs/project.yaml` (`compute.remote_vm`) is
informational only, is a placeholder to configure per site, and is not required to
run any stage.

## (d) GNINA via Docker

GNINA rescoring can run from a local binary or the official `gnina/gnina`
Docker image. Set `gnina.use_docker: true` (Stage 4 / Stage 8 configs) to use:

```
docker run --rm --gpus all -v $PWD:/work gnina/gnina:latest gnina <args>
```

## (e) WSL2 disk space for MD trajectories

The finalist MD run (6 finalists x 3 replicates = 18 x 20 ns) produces ~15 GB of
raw trajectories. Analysis trajectories are water-stripped (~20x smaller), but if
the WSL2 disk (`/`, an ext4 vhdx) runs low, expand it from Windows (plenty free
on C:). This requires shutting WSL down (stops all containers + this SSH session).

From an **elevated PowerShell on Windows**:

```powershell
wsl --shutdown
# find the distro's vhdx (usually under %LOCALAPPDATA%\Packages\...\LocalState\ext4.vhdx)
# grow the max size to e.g. 400 GB:
wsl --manage <DistroName> --set-sparse false           # optional: reclaim sparse
diskpart
  select vdisk file="C:\path\to\ext4.vhdx"
  expand vdisk maximum=409600
  exit
```

Then start WSL and grow the ext4 filesystem online:

```bash
sudo mount -t devtmpfs none /dev 2>/dev/null || true
sudo resize2fs /dev/sdd            # device from `df -T /`
df -h /
```

The MD driver (`scripts/run_finalist_md.sh`) guards against filling the disk: it
aborts before the production phase if free space < EGFR_MD_MIN_FREE_GB (default
4 GB) and resumes cleanly once space is available.
