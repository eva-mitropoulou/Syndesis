"""Virtual-screening enrichment metrics with bootstrap confidence intervals.

Given a scored table (one row per ligand, a real-valued score where higher = more
active-like, and a binary label 1=active/0=decoy), compute the standard LIT-PCBA
metrics: ROC-AUC, enrichment factor at 1%/5% (EF1/EF5), and BEDROC (alpha=80.5).
All are reported with percentile bootstrap 95% CIs by resampling ligands.

Pure/deterministic (seeded bootstrap indices are caller-provided or derived from a
fixed seed) so results are reproducible and unit-testable without any docking.
"""
from __future__ import annotations

from typing import Any

import numpy as np
import pandas as pd


def roc_auc(labels: np.ndarray, scores: np.ndarray) -> float:
    """ROC-AUC via the Mann-Whitney U relation. Ties handled by average rank."""
    labels = np.asarray(labels)
    scores = np.asarray(scores, dtype=float)
    n_pos = int((labels == 1).sum())
    n_neg = int((labels == 0).sum())
    if n_pos == 0 or n_neg == 0:
        return float("nan")
    order = np.argsort(scores, kind="mergesort")
    ranks = np.empty(len(scores), dtype=float)
    sorted_scores = scores[order]
    # average ranks for ties
    i = 0
    r = 1
    while i < len(sorted_scores):
        j = i
        while j < len(sorted_scores) and sorted_scores[j] == sorted_scores[i]:
            j += 1
        avg = (r + (r + (j - i) - 1)) / 2.0
        ranks[order[i:j]] = avg
        r += (j - i)
        i = j
    sum_ranks_pos = ranks[labels == 1].sum()
    auc = (sum_ranks_pos - n_pos * (n_pos + 1) / 2.0) / (n_pos * n_neg)
    return float(auc)


def enrichment_factor(labels: np.ndarray, scores: np.ndarray, fraction: float) -> float:
    """EF at the top `fraction` of the ranked list.
    EF = (actives in top fraction / N_top) / (total actives / N)."""
    labels = np.asarray(labels)
    scores = np.asarray(scores, dtype=float)
    n = len(labels)
    n_top = max(1, int(round(n * fraction)))
    total_active_rate = (labels == 1).mean()
    if total_active_rate == 0:
        return float("nan")
    order = np.argsort(-scores, kind="mergesort")  # descending: best first
    top = labels[order][:n_top]
    top_active_rate = (top == 1).mean()
    return float(top_active_rate / total_active_rate)


def bedroc(labels: np.ndarray, scores: np.ndarray, alpha: float = 80.5) -> float:
    """BEDROC (Truchon & Bayly 2007). alpha=80.5 => ~80% of the score comes from
    the top ~2% of the ranking (the LIT-PCBA convention)."""
    labels = np.asarray(labels)
    scores = np.asarray(scores, dtype=float)
    n = len(labels)
    n_pos = int((labels == 1).sum())
    if n_pos == 0 or n_pos == n:
        return float("nan")
    order = np.argsort(-scores, kind="mergesort")
    ranks = np.where(labels[order] == 1)[0] + 1  # 1-based ranks of actives
    ra = n_pos / n  # fraction of actives
    # RIE = mean_i exp(-alpha * rank_i / N) / expected-value-under-uniform
    sum_exp = np.exp(-alpha * ranks / n).sum()
    rie_denom = (n_pos / n) * (1 - np.exp(-alpha)) / (np.exp(alpha / n) - 1)
    rie = sum_exp / rie_denom
    # BEDROC = RIE * ra*sinh(alpha/2) / (cosh(alpha/2) - cosh(alpha/2 - alpha*ra))
    #          + 1 / (1 - exp(alpha*(1-ra)))
    scale = ra * np.sinh(alpha / 2.0) / (np.cosh(alpha / 2.0) - np.cosh(alpha / 2.0 - alpha * ra))
    bedroc_val = rie * scale + 1.0 / (1.0 - np.exp(alpha * (1.0 - ra)))
    return float(bedroc_val)


def _bootstrap_ci(func, labels: np.ndarray, scores: np.ndarray, n_boot: int, seed: int) -> tuple[float, float]:
    rng = np.random.default_rng(seed)
    n = len(labels)
    vals = []
    for _ in range(n_boot):
        idx = rng.integers(0, n, n)
        v = func(labels[idx], scores[idx])
        if not np.isnan(v):
            vals.append(v)
    if not vals:
        return (float("nan"), float("nan"))
    return (float(np.percentile(vals, 2.5)), float(np.percentile(vals, 97.5)))


def enrichment_report(
    df: pd.DataFrame,
    score_col: str,
    label_col: str = "label",
    *,
    n_boot: int = 2000,
    seed: int = 807,
) -> dict[str, Any]:
    """Full enrichment report for one scoring arm, with bootstrap 95% CIs.
    Rows with a missing score are dropped (and counted) so a failed dock/rescore
    doesn't silently count as a decoy-like 0."""
    sub = df[[score_col, label_col]].copy()
    n_total = len(sub)
    sub = sub.dropna(subset=[score_col])
    n_scored = len(sub)
    labels = sub[label_col].to_numpy()
    scores = sub[score_col].to_numpy(dtype=float)
    out: dict[str, Any] = {
        "arm": score_col,
        "n_total": n_total,
        "n_scored": n_scored,
        "n_actives": int((labels == 1).sum()),
        "n_decoys": int((labels == 0).sum()),
    }
    for name, fn in [
        ("roc_auc", roc_auc),
        ("ef1", lambda l, s: enrichment_factor(l, s, 0.01)),
        ("ef5", lambda l, s: enrichment_factor(l, s, 0.05)),
        ("bedroc", bedroc),
    ]:
        point = fn(labels, scores)
        lo, hi = _bootstrap_ci(fn, labels, scores, n_boot, seed)
        out[name] = round(point, 4) if not np.isnan(point) else None
        out[f"{name}_ci_lo"] = round(lo, 4) if not np.isnan(lo) else None
        out[f"{name}_ci_hi"] = round(hi, 4) if not np.isnan(hi) else None
    return out
