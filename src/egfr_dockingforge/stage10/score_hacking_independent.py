"""De-circularized score-hacking detection with an independent oracle.

The original ``score_hacking`` (metric_definitions.py) and the acceptance gate
``accepted_pre_md`` both key on ``binding_mode_preserved_flag`` and
``delta_candidate_score``, so "accepted" and "score-hacked" are near-complementary
by construction and the reported 0.0 hard rate is close to tautological.

This module decouples them:

* ORACLE (ground truth for "did the pose actually keep the native binding mode"):
  the continuous change in ProLIF IFP-Tanimoto to the native consensus relative to
  the parent (``delta_ifp_tanimoto_to_consensus``). This raw structural-fidelity
  measure does NOT feed the acceptance gate. A molecule whose IFP similarity to the
  validated binding mode drops materially has lost the mode, regardless of any flag.

* DETECTOR (independent of the gate): flags score-hacking when a docking/CNN score
  improves (``delta_gnina_cnnscore > 0`` or ``delta_candidate_score > 0``) while the
  ORACLE says interaction fidelity dropped. Uses only continuous score deltas + the
  oracle, never ``binding_mode_preserved_flag``.

We then report the detector's PRECISION/RECALL against the oracle, so resistance
to score-hacking is measured detector performance, not a definitional artifact.
An adversarial docking-score-only arm
(select purely by score) is the positive control: it should trip the detector.
"""
from __future__ import annotations

from typing import Any

import pandas as pd


# A material loss of native-interaction fidelity (Tanimoto). Below this, changes are
# treated as noise and not called a lost binding mode.
DEFAULT_FIDELITY_DROP = 0.05


def _delta_ifp(row: dict) -> float | None:
    """Change in native-interaction fidelity vs the PARENT (the independent oracle
    signal). This is decoupled from the acceptance gate, which keys on the boolean
    ``binding_mode_preserved_flag`` rather than these continuous deltas.

    Preference order (all measure the same thing — loss of native contacts):
      1. explicit IFP-Tanimoto delta if present;
      2. best_ - parent_ IFP-Tanimoto if both present;
      3. ``delta_key_interaction_recall`` — the analog's key-interaction recall
         minus the parent's; always present in the Stage-10 master and NOT used by
         accepted_pre_md, so it is a valid independent fidelity oracle.
    """
    d = row.get("delta_ifp_tanimoto_to_consensus")
    if d is not None:
        return float(d)
    best = row.get("best_ifp_tanimoto_to_consensus")
    parent = row.get("parent_ifp_tanimoto_to_consensus")
    if best is not None and parent is not None:
        return float(best) - float(parent)
    dk = row.get("delta_key_interaction_recall")
    if dk is not None:
        return float(dk)
    return None


def oracle_lost_binding_mode(row: dict, fidelity_drop: float = DEFAULT_FIDELITY_DROP) -> bool | None:
    """Independent ground truth: did the analog lose the native binding mode?
    True if IFP-Tanimoto-to-consensus dropped by more than `fidelity_drop`, False if
    it held, None if the oracle signal is unavailable."""
    d = _delta_ifp(row)
    if d is None:
        return None
    return d < -abs(fidelity_drop)


def score_improved(row: dict) -> bool:
    """A docking/CNN score got better (the thing a hacker optimizes)."""
    return bool(row.get("delta_gnina_cnnscore", 0.0) > 0 or row.get("delta_candidate_score", 0.0) > 0)


def score_hacking_independent(row: dict, fidelity_drop: float = DEFAULT_FIDELITY_DROP) -> bool | None:
    """Score-hacking flag using ONLY continuous score deltas + the independent
    oracle: score improved while the oracle says interaction fidelity dropped.
    None if the oracle signal is unavailable for this row."""
    lost = oracle_lost_binding_mode(row, fidelity_drop)
    if lost is None:
        return None
    return bool(score_improved(row) and lost)


def evaluate_detector(master: pd.DataFrame, fidelity_drop: float = DEFAULT_FIDELITY_DROP) -> dict[str, Any]:
    """Precision/recall of the OLD gate-coupled score_hacking flag against the new
    independent oracle, plus per-strategy independent hacking rates. Rows without
    an oracle signal are excluded and counted.

    The oracle positive class = "true score-hacker" := score_improved AND oracle
    lost_binding_mode. Precision/recall are computed for the ORIGINAL
    ``score_hacking_flag`` (the gate-coupled detector) treated as the prediction,
    so we can state how well the reported metric tracks the independent truth."""
    from egfr_dockingforge.stage10.metric_definitions import score_hacking as gate_coupled

    rows = master.to_dict("records")
    tp = fp = fn = tn = 0
    n_oracle = 0
    ind_flags = []
    for r in rows:
        truth = score_hacking_independent(r, fidelity_drop)
        if truth is None:
            ind_flags.append(None)
            continue
        n_oracle += 1
        ind_flags.append(truth)
        pred = bool(r.get("score_hacking_flag")) if "score_hacking_flag" in r else bool(gate_coupled(r))
        if pred and truth:
            tp += 1
        elif pred and not truth:
            fp += 1
        elif not pred and truth:
            fn += 1
        else:
            tn += 1
    precision = tp / (tp + fp) if (tp + fp) else None
    recall = tp / (tp + fn) if (tp + fn) else None
    out: dict[str, Any] = {
        "n_total": len(rows),
        "n_with_oracle": n_oracle,
        "n_independent_hackers": int(sum(1 for f in ind_flags if f)),
        "gate_detector_precision_vs_oracle": None if precision is None else round(precision, 3),
        "gate_detector_recall_vs_oracle": None if recall is None else round(recall, 3),
        "confusion": {"tp": tp, "fp": fp, "fn": fn, "tn": tn},
        "fidelity_drop_threshold": fidelity_drop,
    }
    return out


def threshold_sweep(master: pd.DataFrame, thresholds: tuple[float, ...] = (0.01, 0.02, 0.05)) -> pd.DataFrame:
    """Detector behaviour vs the independent oracle across fidelity-drop thresholds.
    Reporting the sweep (not one number) is the honest presentation: it shows how
    many 'score improved but lost N native-interaction bits' cases exist at each
    severity, and how well the original gate-coupled flag tracks them."""
    rows = []
    for t in thresholds:
        ev = evaluate_detector(master, fidelity_drop=t)
        rows.append({
            "fidelity_drop_threshold": t,
            "n_independent_hackers": ev["n_independent_hackers"],
            "gate_precision_vs_oracle": ev["gate_detector_precision_vs_oracle"],
            "gate_recall_vs_oracle": ev["gate_detector_recall_vs_oracle"],
            "tp": ev["confusion"]["tp"], "fp": ev["confusion"]["fp"],
            "fn": ev["confusion"]["fn"], "tn": ev["confusion"]["tn"],
        })
    return pd.DataFrame(rows)


def run_independent_score_hacking(master: pd.DataFrame, out_dir) -> dict[str, Any]:
    """Full de-circularized analysis: threshold sweep + per-strategy independent
    rates, written to disk. Returns a summary dict."""
    from pathlib import Path
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    sweep = threshold_sweep(master)
    per_strat = per_strategy_independent_rate(master)
    sweep.to_csv(out_dir / "score_hacking_independent_sweep.csv", index=False)
    sweep.to_parquet(out_dir / "score_hacking_independent_sweep.parquet", index=False)
    per_strat.to_csv(out_dir / "score_hacking_independent_by_strategy.csv", index=False)
    per_strat.to_parquet(out_dir / "score_hacking_independent_by_strategy.parquet", index=False)
    return {"sweep": sweep.to_dict("records"), "by_strategy": per_strat.to_dict("records")}


def per_strategy_independent_rate(master: pd.DataFrame, fidelity_drop: float = DEFAULT_FIDELITY_DROP) -> pd.DataFrame:
    """Independent score-hacking rate per strategy (fraction of oracle-evaluable
    analogs flagged by the independent detector)."""
    m = master.copy()
    m["_ind"] = m.apply(lambda r: score_hacking_independent(r.to_dict(), fidelity_drop), axis=1)
    rows = []
    for strat, g in m.groupby("strategy_name"):
        evaluable = g[g["_ind"].notna()]
        rows.append({
            "strategy_name": strat,
            "n_analogs": len(g),
            "n_oracle_evaluable": len(evaluable),
            "independent_score_hacking_rate": (float(evaluable["_ind"].mean()) if len(evaluable) else None),
        })
    return pd.DataFrame(rows)
