import pandas as pd

from egfr_dockingforge.stage11.gromacs_runner import make_gromacs_runs


def test_gromacs_command_rows_record_blocked_status(tmp_path):
    systems = pd.DataFrame([{"md_system_id":"s1","md_candidate_id":"c1","ionized_structure_file":"","complex_file":"","topology_file":"","build_status":"blocked_missing_ligand_parameters"}])
    config = {"forcefield":{"gromacs_executable":"/usr/local/gromacs/bin/gmx"},"md":{"gpu":"auto"}}
    out = make_gromacs_runs(systems, config, {"processed": tmp_path, "md_root": tmp_path})
    assert not out.empty
    assert out["command_line_grompp"].str.contains("grompp").all()
    assert out["run_status"].eq("blocked_missing_system_build").all()
