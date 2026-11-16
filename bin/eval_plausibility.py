#!/usr/bin/env python3
"""Thin CLI wrapper around soillifeatlas.evaluation.plausibility_score.

Scores an observed kingdom-composition against expected grassland-soil
ranges from the literature. Collapses the rich dict returned by the
framework into a single-row TSV per method invocation — so the pipeline
can concatenate multiple method outputs into one flat table for the CI
gate.

Inputs:
  --composition-kingdom  Parquet in long format (sample_id, kingdom,
                         proportion_pct). Aggregated to mean across samples
                         per kingdom before scoring — plausibility is a
                         cohort-level metric.
  --expected-ref         CSV with columns kingdom, expected_min_pct,
                         expected_max_pct, expected_midpoint_pct
  --method               Method tag emitted into the output (e.g.
                         "fc_weighted_bc"). Free-form; not validated.
  --output               Output TSV, one row per invocation

Output columns (contract):
  method, bc_vs_expected, in_range_fraction, n_kingdoms_scored,
  inflation_score, deviation_<kingdom>*

Note: `bc_vs_expected` is the emitted alias for the upstream framework's
`bc_distance_from_midpoint` key, so the CI assertion
    plaus[plaus.method == "fc_weighted_bc"]["bc_vs_expected"] <= 0.15
lines up without callers having to know the internal naming.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

from soillifeatlas.evaluation import plausibility_score


def _flatten_score(method: str, score: dict) -> dict:
    """Flatten the framework's nested-dict score into a single flat row.

    The upstream dict carries:
      - bc_distance_from_midpoint           -> bc_vs_expected (renamed for CI)
      - in_range_fraction                   -> in_range_fraction
      - n_kingdoms_scored                   -> n_kingdoms_scored
      - inflation_score_sum                 -> inflation_score (renamed)
      - deviation_per_kingdom: {k: float}   -> deviation_<k> columns
      - (observed_pct, expected_midpoint_pct, in_range_per_kingdom,
         deviation_per_kingdom): per-kingdom dicts we mostly drop, keeping
         only the deviations since they're what's informative when the BC
         is too high.
    """
    row: dict = {
        "method": method,
        "bc_vs_expected": float(score["bc_distance_from_midpoint"]),
        "in_range_fraction": float(score["in_range_fraction"]),
        "n_kingdoms_scored": int(score["n_kingdoms_scored"]),
        "inflation_score": float(score["inflation_score_sum"])
        if pd.notna(score["inflation_score_sum"]) else float("nan"),
    }
    for kingdom, dev in score.get("deviation_per_kingdom", {}).items():
        row[f"deviation_{kingdom}"] = float(dev)
    return row


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--composition-kingdom", required=True, type=Path,
                   help="Long-format kingdom composition parquet "
                        "(sample_id, kingdom, proportion_pct)")
    p.add_argument("--expected-ref", required=True, type=Path,
                   help="Expected kingdom composition CSV "
                        "(kingdom, expected_min_pct, expected_max_pct, expected_midpoint_pct)")
    p.add_argument("--method", required=True,
                   help="Method tag to emit in the output (e.g. fc_weighted_bc)")
    p.add_argument("--output", required=True, type=Path,
                   help="Output TSV (one row per invocation)")
    args = p.parse_args(argv)

    # --- Load inputs ---------------------------------------------------------
    comp = pd.read_parquet(args.composition_kingdom)
    ref = pd.read_csv(args.expected_ref)

    required_cols = {"sample_id", "kingdom", "proportion_pct"}
    if not required_cols.issubset(comp.columns):
        sys.exit(
            f"composition_kingdom is missing required columns; "
            f"got {list(comp.columns)}, need {sorted(required_cols)}"
        )

    # --- Aggregate to mean-across-samples series (one value per kingdom) ----
    observed_kingdom_pct = (
        comp.groupby("kingdom")["proportion_pct"]
            .mean()
            .astype(float)
    )

    # --- Score ---------------------------------------------------------------
    score = plausibility_score(
        observed_kingdom_pct=observed_kingdom_pct,
        expected_ref=ref,
    )

    # --- Flatten and emit ----------------------------------------------------
    row = _flatten_score(args.method, score)
    out_df = pd.DataFrame([row])
    args.output.parent.mkdir(parents=True, exist_ok=True)
    out_df.to_csv(args.output, sep="\t", index=False)
    return 0


if __name__ == "__main__":
    sys.exit(main())
