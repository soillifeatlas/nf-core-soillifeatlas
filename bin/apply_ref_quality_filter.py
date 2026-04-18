#!/usr/bin/env python3
"""Apply Layer 5 reference-quality filter to a SIMPER fingerprint atlas.

Restricts Archaea reference rows to the ArchLips-validated feature set, as
described by CorrectionConfig.restrict_archaea_to_archlips. Non-Archaea rows
(Bacteria, Fungi, ...) are passed through unchanged.

CorrectionConfig.restrict_archaea_to_archlips is defined on the dataclass but
is NOT implemented inside corrections.correct_intensity (which only applies
biomass / IS / RIE / external calibration). The filter therefore runs as a
direct atlas-level mask in this wrapper.

Inputs:
  --simper-atlas        SIMPER fingerprint atlas parquet (requires columns
                        `feature_id`, `kingdom`)
  --archlips-validated  CSV with column `feature_id` listing ArchLips-validated
                        features
  --archaea-kingdom     Value used in the atlas `kingdom` column for archaea
                        rows (default: "Archaea")
  --output              Filtered atlas parquet
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd


def apply_ref_quality_filter(
    atlas: pd.DataFrame,
    validated_feature_ids: set[str],
    archaea_kingdom: str = "Archaea",
) -> pd.DataFrame:
    """Drop atlas rows where `kingdom == archaea_kingdom` and `feature_id` is
    not in `validated_feature_ids`. Non-archaea rows are preserved verbatim.
    """
    required = {"feature_id", "kingdom"}
    missing = required - set(atlas.columns)
    if missing:
        raise ValueError(
            f"SIMPER atlas must have columns {sorted(required)}, missing: {sorted(missing)}"
        )
    is_archaea = atlas["kingdom"] == archaea_kingdom
    keep = ~is_archaea | atlas["feature_id"].isin(validated_feature_ids)
    return atlas.loc[keep].reset_index(drop=True)


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--simper-atlas", required=True, type=Path,
                   help="SIMPER fingerprint atlas parquet")
    p.add_argument("--archlips-validated", required=True, type=Path,
                   help="ArchLips-validated feature CSV (column: feature_id)")
    p.add_argument("--archaea-kingdom", default="Archaea",
                   help='Kingdom value marking archaea rows (default: "Archaea")')
    p.add_argument("--output", required=True, type=Path,
                   help="Filtered SIMPER atlas parquet")
    args = p.parse_args(argv)

    atlas = pd.read_parquet(args.simper_atlas)
    archlips = pd.read_csv(args.archlips_validated)
    if "feature_id" not in archlips.columns:
        raise ValueError(
            f"archlips-validated CSV must have a `feature_id` column "
            f"(got columns: {list(archlips.columns)})"
        )
    validated = set(archlips["feature_id"])

    filtered = apply_ref_quality_filter(atlas, validated, archaea_kingdom=args.archaea_kingdom)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    filtered.to_parquet(args.output)
    return 0


if __name__ == "__main__":
    sys.exit(main())
