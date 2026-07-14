from pathlib import Path

from egfr_dockingforge.stage11.system_builder import _ionized_system_counts


def test_ionized_system_counts_are_read_from_topology(tmp_path: Path) -> None:
    topology = tmp_path / "topol.top"
    topology.write_text(
        """
[ system ]
Protein in water

[ molecules ]
Protein_chain_A 1
UNL 1
SOL 32285
NA 103
CL 95
"""
    )

    assert _ionized_system_counts(topology) == {
        "num_waters": 32285,
        "num_na": 103,
        "num_cl": 95,
        "net_charge_before_ions": -8,
        "final_charge": 0,
    }
