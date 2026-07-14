from __future__ import annotations

import numpy as np
import pandas as pd
from sklearn.calibration import CalibratedClassifierCV, calibration_curve
from sklearn.frozen import FrozenEstimator
from sklearn.metrics import brier_score_loss

from syndesis.common.io import ensure_dir, write_table


def expected_calibration_error(y_true: np.ndarray, y_prob: np.ndarray, bins: int = 10) -> float:
    edges = np.linspace(0, 1, bins + 1)
    ece = 0.0
    for lo, hi in zip(edges[:-1], edges[1:], strict=False):
        mask = (y_prob >= lo) & (y_prob < hi if hi < 1 else y_prob <= hi)
        if not np.any(mask):
            continue
        ece += mask.mean() * abs(float(y_true[mask].mean()) - float(y_prob[mask].mean()))
    return float(ece)


def calibrate_classifier(model: object, x_valid, y_valid, model_id: str, paths: dict) -> tuple[object, pd.DataFrame]:
    method = "isotonic" if len(y_valid) >= 80 and len(set(y_valid)) == 2 else "sigmoid"
    before = model.predict_proba(x_valid)[:, 1]
    calibrated = CalibratedClassifierCV(FrozenEstimator(model), method=method)
    calibrated.fit(x_valid, y_valid)
    after = calibrated.predict_proba(x_valid)[:, 1]
    curve = pd.DataFrame()
    if len(set(y_valid)) == 2:
        prob_true, prob_pred = calibration_curve(y_valid, after, n_bins=10, strategy="uniform")
        curve = pd.DataFrame({"prob_true": prob_true, "prob_pred": prob_pred})
    curve_path = paths["processed"] / f"reliability_curve_{model_id}.csv"
    ensure_dir(curve_path.parent)
    curve.to_csv(curve_path, index=False)
    metrics = pd.DataFrame(
        [
            {
                "model_id": model_id,
                "calibration_method": method,
                "split_name": "primary",
                "brier_score_before": brier_score_loss(y_valid, before),
                "brier_score_after": brier_score_loss(y_valid, after),
                "ece_before": expected_calibration_error(np.asarray(y_valid), before),
                "ece_after": expected_calibration_error(np.asarray(y_valid), after),
                "reliability_curve_path": str(curve_path),
                "calibration_status": "complete",
            }
        ]
    )
    write_table(paths["processed"] / "calibration_metrics.parquet", metrics)
    write_table(paths["processed"] / "calibration_metrics.csv", metrics)
    return calibrated, metrics
