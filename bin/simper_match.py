#!/usr/bin/env python3
"""MS2 cosine matching of soil spectra vs atlas spectra with SIMPER filter.

This is the critical bridge from soil MS2 data to the atlas SIMPER
decomposition reference. Per Samrat 2025 methods:

  - Modified cosine similarity (matchms ``ModifiedCosineGreedy``)
  - Precursor tolerance: 10 ppm
  - Fragment tolerance: 0.02 Da
  - ``top_n = 50`` peaks kept per spectrum
  - Cosine threshold ``>= 0.7``
  - Minimum 4 matched peaks

Given (a) soil MS2 MGF, (b) atlas MS2 MGF, and (c) the SIMPER fingerprint
atlas parquet, emit a long-format parquet of verified matches joined to
the SIMPER fingerprint (phylum, kingdom, direction, fold_change,
simper_rank). One row per (soil_feature_id, atlas_feature_id, phylum)
passing filters — a single soil feature may match multiple atlas features
across multiple phyla; downstream decomposition handles the many-to-many.

matchms API notes (0.32.0):
  - The "modified cosine" operator was renamed. Use
    ``matchms.similarity.ModifiedCosineGreedy`` (greedy peak assignment;
    this is the standard, fast variant used in the paper). The exact
    Hungarian-assignment variant ``ModifiedCosineHungarian`` also exists
    but is slower.
  - ``mc.pair(ref, query)`` returns a structured ``numpy.ndarray`` with
    dtype ``[('score', '<f8'), ('matches', '<i8')]``. Access via
    ``res['score']`` and ``res['matches']``.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from matchms.importing import load_from_mgf
from matchms.filtering import (
    default_filters,
    normalize_intensities,
    reduce_to_number_of_peaks,
)
from matchms.similarity import ModifiedCosineGreedy


# ---------------------------------------------------------------------------
# Output schema — exported so downstream tooling can align empty DataFrames
# ---------------------------------------------------------------------------
OUTPUT_COLUMNS = [
    "soil_feature_id",
    "atlas_feature_id",
    "phylum",
    "kingdom",
    "direction",
    "fold_change",
    "simper_rank",
    "cosine",
    "n_matched_peaks",
    "precursor_ppm_diff",
]

SIMPER_KEEP_COLS = [
    "feature_id", "phylum", "kingdom", "direction",
    "fold_change", "simper_rank",
]


def _preprocess(spec, top_n: int):
    """Apply the matchms preprocessing chain required for ModifiedCosine.

    matchms strongly recommends (a) normalizing intensities and (b) capping
    the peak count before cosine scoring. ``default_filters`` also fixes
    metadata keys (e.g. PEPMASS -> precursor_mz) that vary across MGF
    producers.
    """
    if spec is None:
        return None
    spec = default_filters(spec)
    if spec is None:
        return None
    spec = normalize_intensities(spec)
    if spec is None:
        return None
    spec = reduce_to_number_of_peaks(spec, n_max=top_n)
    return spec


def _extract_feature_id(spec) -> str | None:
    """Pull a stable feature_id from a matchms Spectrum.

    Prefer ``FEATURE_ID`` (custom MGF field we emit), fall back to
    ``SCANS`` (standard MGF field), then ``TITLE``.
    """
    return (
        spec.get("feature_id")
        or spec.get("scans")
        or spec.get("title")
    )


def _get_precursor(spec) -> float | None:
    """Get precursor m/z as float, or None if missing."""
    prec = spec.get("precursor_mz")
    if prec is None:
        prec = spec.get("pepmass")
        # pepmass can be a tuple (mz, intensity); take the first
        if isinstance(prec, (tuple, list)) and len(prec) > 0:
            prec = prec[0]
    if prec is None:
        return None
    try:
        return float(prec)
    except (TypeError, ValueError):
        return None


def _empty_output() -> pd.DataFrame:
    """Return an empty DataFrame with the full 10-column schema."""
    return pd.DataFrame({c: pd.Series(dtype=object) for c in OUTPUT_COLUMNS})


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--soil-mgf", required=True, type=Path,
                   help="Soil MS2 spectra in MGF format")
    p.add_argument("--atlas-mgf", required=True, type=Path,
                   help="Atlas (consensus) MS2 spectra in MGF format")
    p.add_argument("--simper-atlas", required=True, type=Path,
                   help="SIMPER fingerprint atlas parquet (feature_id, "
                        "phylum, kingdom, direction, fold_change, simper_rank)")
    p.add_argument("--precursor-ppm", type=float, default=10.0,
                   help="Precursor m/z tolerance in ppm (default 10)")
    p.add_argument("--fragment-tol", type=float, default=0.02,
                   help="Fragment m/z tolerance in Da (default 0.02)")
    p.add_argument("--top-n", type=int, default=50,
                   help="Keep top-N peaks per spectrum (default 50)")
    p.add_argument("--min-cos", type=float, default=0.7,
                   help="Minimum modified-cosine score to accept (default 0.7)")
    p.add_argument("--min-matched-peaks", type=int, default=4,
                   help="Minimum number of matched peaks (default 4)")
    p.add_argument("--output-matches", required=True, type=Path,
                   help="Output parquet of verified matches (long format)")
    args = p.parse_args(argv)

    # --- Load + preprocess ------------------------------------------------
    soil_spectra = [_preprocess(s, args.top_n)
                    for s in load_from_mgf(str(args.soil_mgf))]
    atlas_spectra = [_preprocess(s, args.top_n)
                     for s in load_from_mgf(str(args.atlas_mgf))]
    soil_spectra = [s for s in soil_spectra
                    if s is not None and len(s.peaks.mz) > 0]
    atlas_spectra = [a for a in atlas_spectra
                     if a is not None and len(a.peaks.mz) > 0]

    # Pre-extract atlas feature_ids + precursors to avoid re-parsing per
    # soil spectrum. This is the inner-loop optimization that matters when
    # M (atlas spectra) is in the hundreds of thousands.
    atlas_index: list[tuple[object, float, str]] = []
    for a in atlas_spectra:
        a_prec = _get_precursor(a)
        a_fid = _extract_feature_id(a)
        if a_prec is None or a_fid is None:
            continue
        atlas_index.append((a, a_prec, a_fid))

    mc = ModifiedCosineGreedy(tolerance=args.fragment_tol)

    rows: list[dict] = []
    for s_spec in soil_spectra:
        s_fid = _extract_feature_id(s_spec)
        if s_fid is None:
            continue
        s_prec = _get_precursor(s_spec)
        if s_prec is None:
            continue

        for a_spec, a_prec, a_fid in atlas_index:
            # Precursor tolerance check (cheap; gate before cosine math)
            if a_prec <= 0:
                continue
            prec_ppm_diff = abs(s_prec - a_prec) / a_prec * 1e6
            if prec_ppm_diff > args.precursor_ppm:
                continue

            score = mc.pair(s_spec, a_spec)
            cos = float(score["score"])
            n_matched = int(score["matches"])
            if cos < args.min_cos or n_matched < args.min_matched_peaks:
                continue

            rows.append({
                "soil_feature_id": s_fid,
                "atlas_feature_id": a_fid,
                "cosine": cos,
                "n_matched_peaks": n_matched,
                "precursor_ppm_diff": prec_ppm_diff,
            })

    matches_df = pd.DataFrame(rows)

    # --- Join against SIMPER atlas ---------------------------------------
    simper = pd.read_parquet(args.simper_atlas)
    simper_slim = (
        simper[[c for c in SIMPER_KEEP_COLS if c in simper.columns]]
        .drop_duplicates()
    )

    if len(matches_df) > 0:
        joined = matches_df.merge(
            simper_slim,
            left_on="atlas_feature_id",
            right_on="feature_id",
            how="inner",
        ).drop(columns=["feature_id"])
        # Re-order columns to the published schema (stable contract).
        joined = joined[[c for c in OUTPUT_COLUMNS if c in joined.columns]]
    else:
        joined = _empty_output()

    args.output_matches.parent.mkdir(parents=True, exist_ok=True)
    joined.to_parquet(args.output_matches)

    # Terse pipeline-log summary to stderr
    print(
        f"simper_match: soil_spectra={len(soil_spectra)} "
        f"atlas_spectra={len(atlas_spectra)} "
        f"raw_matches={len(matches_df)} verified_rows={len(joined)}",
        file=sys.stderr,
    )
    return 0


if __name__ == "__main__":
    sys.exit(main() or 0)
