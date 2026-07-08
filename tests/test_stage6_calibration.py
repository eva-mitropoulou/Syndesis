from __future__ import annotations

import numpy as np

from egfr_dockingforge.stage6.calibration import expected_calibration_error


def test_expected_calibration_error_in_unit_interval() -> None:
    value = expected_calibration_error(np.array([0, 1, 1, 0]), np.array([0.1, 0.8, 0.7, 0.2]))
    assert 0 <= value <= 1
