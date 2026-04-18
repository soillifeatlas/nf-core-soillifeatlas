#!/usr/bin/env python3
"""Thin CLI wrapper around soillifeatlas.corrections.apply_RIE_correction.

Layer 3 of the DECOMPOSE_APPLY subworkflow: per-feature scaling by the class +
adduct relative ionization efficiency (RIE). A floor/ceiling caps extreme
amplification from very-low or very-high RIE values.

Inputs:
  --intensity    Parquet (features x samples, feature_id as index)
  --rie-table    CSV with columns: class, adduct, RIE_LPE
  --annotation   CSV with columns: feature_id, class, adduct
  --rie-floor    Min RIE allowed (e.g. 0.05 -> max 20x amplification; 0.20 -> max 5x)
  --rie-ceiling  Max RIE allowed (default 100.0, matches CorrectionConfig)
  --fallback-rie RIE for features missing from the annotation/RIE table (default 1.0)
  --output       Output RIE-corrected parquet
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

from soillifeatlas.corrections import apply_RIE_correction


def _align_annotation(intensity: pd.DataFrame, annotation: pd.DataFrame) -> tuple[pd.Series, pd.Series]:
    """Return (feature_class, feature_adduct) indexed to match `intensity.index`.

    Rows in `intensity` without an annotation row get NaN for class/adduct, which
    the underlying apply_RIE_correction treats as fallback_RIE.
    """
    required = {"feature_id", "class", "adduct"}
    missing = required - set(annotation.columns)
    if missing:
        raise ValueError(
            f"annotation CSV must have columns {sorted(required)}, missing: {sorted(missing)}"
        )
    ann = annotation.set_index("feature_id").reindex(intensity.index)
    return ann["class"], ann["adduct"]


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--intensity", required=True, type=Path,
                   help="Input intensity parquet (features x samples)")
    p.add_argument("--rie-table", required=True, type=Path,
                   help="RIE table CSV (columns: class, adduct, RIE_LPE)")
    p.add_argument("--annotation", required=True, type=Path,
                   help="Feature annotation CSV (columns: feature_id, class, adduct)")
    p.add_argument("--rie-floor", type=float, required=True,
                   help="Minimum RIE; caps amplification at 1/floor")
    p.add_argument("--rie-ceiling", type=float, default=100.0,
                   help="Maximum RIE; caps suppression at 1/ceiling (default: 100.0)")
    p.add_argument("--fallback-rie", type=float, default=1.0,
                   help="RIE for features without class/adduct annotation (default: 1.0)")
    p.add_argument("--output", required=True, type=Path,
                   help="Output RIE-corrected parquet")
    args = p.parse_args(argv)

    intensity = pd.read_parquet(args.intensity)
    rie_table = pd.read_csv(args.rie_table)
    annotation = pd.read_csv(args.annotation)

    feature_class, feature_adduct = _align_annotation(intensity, annotation)

    corrected, _rie_vec = apply_RIE_correction(
        intensity=intensity,
        feature_class=feature_class,
        feature_adduct=feature_adduct,
        rie_lookup=rie_table,
        fallback_RIE=args.fallback_rie,
        rie_floor=args.rie_floor,
        rie_ceiling=args.rie_ceiling,
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    corrected.to_parquet(args.output)
    return 0


if __name__ == "__main__":
    sys.exit(main())
