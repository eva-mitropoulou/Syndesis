# Stage 4 GNINA Usage

GNINA is used in score-only mode:

```bash
gnina --score_only -r receptor.pdb -l pose.pdbqt
```

The configured local executable is a Docker wrapper for `gnina/gnina:latest`. The current image reports GNINA v1.3.3. The wrapper mounts the repository root as `/work`, so Stage 4 passes repository-relative input paths whenever possible.

Default policy:

- do not redock
- do not minimize
- do not use covalent mode
- record command line, GNINA version, model, runtime, exit code, warnings, and parse status
- keep missing optional score fields as null with explicit warnings

GNINA scores are treated as rescoring features for Stage 5/6 integration, not final biological evidence.
