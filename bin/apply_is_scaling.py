#!/usr/bin/env python3
"""Thin CLI wrapper around soillifeatlas.corrections.apply_IS_normalization.

Layer 2 of the DECOMPOSE_APPLY subworkflow: per-sample scaling by the internal
standard (IS) intensity. The IS is identified by a compound name that maps to
one feature_id via the IS feature table.

Inputs:
  --intensity     Parquet, features x samples (feature_id as index)
  --is-features   CSV with columns: feature_id, compound, adduct, spiked_pmol
  --is-reference  IS compound name (e.g., LPE_18d7) — must match a row in
                  is_features.csv
  --is-spiked-pmol  Spiked IS amount in pmol per sample
  --output        Output parquet (features x samples, IS-normalized)
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

from soillifeatlas.corrections import apply_IS_normalization


def _resolve_is_row(intensity: pd.DataFrame, is_features: pd.DataFrame, reference_compound: str) -> pd.Series:
    """Return the intensity-matrix row for the IS feature identified by `reference_compound`.

    Raises a clear error if the reference compound is not in is_features or the
    mapped feature_id is not in intensity.
    """
    matches = is_features.loc[is_features["compound"] == reference_compound]
    if matches.empty:
        raise ValueError(
            f"IS reference compound {reference_compound!r} not found in is_features "
            f"(available: {sorted(is_features['compound'].unique())})"
        )
    if len(matches) > 1:
        raise ValueError(
            f"IS reference compound {reference_compound!r} maps to multiple feature_ids: "
            f"{matches['feature_id'].tolist()}"
        )
    is_feature_id = matches.iloc[0]["feature_id"]
    if is_feature_id not in intensity.index:
        raise ValueError(
            f"IS feature_id {is_feature_id!r} (compound={reference_compound!r}) "
            f"is not present in the intensity matrix index"
        )
    return intensity.loc[is_feature_id]


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--intensity", required=True, type=Path,
                   help="Input intensity parquet (features x samples)")
    p.add_argument("--is-features", required=True, type=Path,
                   help="IS feature CSV (columns: feature_id, compound, adduct, spiked_pmol)")
    p.add_argument("--is-spiked-pmol", type=float, required=True,
                   help="Spiked IS amount in pmol per sample (e.g., 100)")
    p.add_argument("--is-reference", required=True,
                   help="IS reference compound name (e.g., LPE_18d7)")
    p.add_argument("--output", required=True, type=Path,
                   help="Output IS-normalized parquet (features x samples)")
    args = p.parse_args(argv)

    intensity = pd.read_parquet(args.intensity)
    is_features = pd.read_csv(args.is_features)

    IS_intensity_per_sample = _resolve_is_row(intensity, is_features, args.is_reference)

    corrected = apply_IS_normalization(
        intensity=intensity,
        IS_intensity_per_sample=IS_intensity_per_sample,
        IS_spiked_pmol=args.is_spiked_pmol,
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    corrected.to_parquet(args.output)
    return 0


if __name__ == "__main__":
    sys.exit(main())
