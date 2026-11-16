"""
Evaluation metrics for quantification pipeline iterations.

Primary metric: biological plausibility of kingdom composition vs grassland
soil literature expectations.

Secondary metric: POS <-> NEG agreement (lower divergence = better).
"""
from __future__ import annotations

from pathlib import Path
import numpy as np
import pandas as pd


def bray_curtis(x: np.ndarray, y: np.ndarray) -> float:
    """Bray-Curtis dissimilarity between two composition vectors (0 = identical, 1 = maximally different)."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    return float(np.sum(np.abs(x - y)) / (np.sum(x) + np.sum(y) + 1e-12))


def plausibility_score(
    observed_kingdom_pct: pd.Series,
    expected_ref: pd.DataFrame,
    exclude_kingdoms: tuple[str, ...] = ("Viruses", "Mixed"),
) -> dict:
    """Score observed kingdom composition vs expected grassland soil ranges.

    Returns dict with:
      - bc_distance_from_midpoint: Bray-Curtis vs expected midpoint (lower = better)
      - in_range_score: fraction of kingdoms with observed in [min, max] (higher = better)
      - deviation_per_kingdom: pd.Series with positive/negative deviation from midpoint
      - inflation_score: absolute deviation on inflated kingdoms (Animalia, Archaea, Protozoa)
                         lower = less inflated = better
    """
    # Filter to common kingdoms
    expected = expected_ref.set_index("kingdom")
    common = [k for k in expected.index if k in observed_kingdom_pct.index and k not in exclude_kingdoms]

    obs = observed_kingdom_pct[common].astype(float)
    # Renormalize to 100% after excluding
    obs = obs / obs.sum() * 100

    mid = expected.loc[common, "expected_midpoint_pct"].astype(float)
    mid = mid / mid.sum() * 100  # also renormalize expected so both sum to 100

    lo = expected.loc[common, "expected_min_pct"].astype(float)
    hi = expected.loc[common, "expected_max_pct"].astype(float)

    bc = bray_curtis(obs.values, mid.values)

    in_range = ((obs >= lo) & (obs <= hi)).astype(int)

    deviation = obs - mid

    # Inflation score: abs deviation for kingdoms currently inflated in raw data
    inflated_kingdoms = ["Animalia", "Archaea", "Protozoa"]
    inflation_targets = [k for k in inflated_kingdoms if k in common]
    # For inflation metric, positive deviation is bad
    inflation = float(deviation[inflation_targets].clip(lower=0).sum()) if inflation_targets else np.nan

    return {
        "bc_distance_from_midpoint": bc,
        "in_range_fraction": float(in_range.mean()),
        "in_range_per_kingdom": in_range.to_dict(),
        "deviation_per_kingdom": deviation.to_dict(),
        "observed_pct": obs.to_dict(),
        "expected_midpoint_pct": mid.to_dict(),
        "inflation_score_sum": inflation,
        "n_kingdoms_scored": len(common),
    }


def pos_neg_agreement(
    pos_composition: pd.Series,
    neg_composition: pd.Series,
    exclude_kingdoms: tuple[str, ...] = ("Viruses", "Mixed"),
) -> dict:
    """Quantify POS <-> NEG kingdom % agreement. Lower BC = better."""
    common = [k for k in pos_composition.index if k in neg_composition.index and k not in exclude_kingdoms]
    p = pos_composition[common].astype(float)
    n = neg_composition[common].astype(float)
    # Renormalize
    p = p / p.sum() * 100
    n = n / n.sum() * 100
    bc = bray_curtis(p.values, n.values)
    max_diff = float((p - n).abs().max())
    return {
        "pos_neg_bc_distance": bc,
        "max_abs_diff_pct": max_diff,
        "per_kingdom_abs_diff": (p - n).abs().to_dict(),
    }


def summarize_iteration(
    iter_name: str,
    observed_pos: pd.Series | None,
    observed_neg: pd.Series | None,
    expected_ref: pd.DataFrame,
) -> dict:
    """One-stop summary of an iteration's quality."""
    out = {"iteration": iter_name}
    if observed_pos is not None:
        out["POS"] = plausibility_score(observed_pos, expected_ref)
    if observed_neg is not None:
        out["NEG"] = plausibility_score(observed_neg, expected_ref)
    if observed_pos is not None and observed_neg is not None:
        out["agreement"] = pos_neg_agreement(observed_pos, observed_neg)
    return out
