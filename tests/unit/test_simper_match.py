"""Unit test for bin/simper_match.py.

TDD: written before the wrapper. Exercises the CLI end-to-end via subprocess
so we verify both arg parsing and the underlying matchms ModifiedCosine call.

This wrapper is the critical bridge from soil features to the decomposition
reference: MS2 cosine matching of soil spectra vs atlas spectra, then filtering
the SIMPER fingerprint atlas to the verified feature_ids.
"""
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from matchms import Spectrum
from matchms.exporting import save_as_mgf


REPO_ROOT = Path(__file__).resolve().parents[2]
BIN = REPO_ROOT / "bin" / "simper_match.py"


def _make_spectrum(feature_id, precursor_mz, peak_mz_list, peak_int_list):
    return Spectrum(
        mz=np.array(peak_mz_list, dtype=float),
        intensities=np.array(peak_int_list, dtype=float),
        metadata={
            "precursor_mz": precursor_mz,
            "feature_id": feature_id,
            "scans": feature_id,
        },
    )


def test_simper_match_end_to_end(tmp_path):
    # Atlas: 2 spectra at known precursors
    atlas = [
        _make_spectrum(
            "atlas_feat_000", 500.123,
            [100.05, 200.10, 300.15, 400.20], [1.0, 0.8, 0.6, 0.4],
        ),
        _make_spectrum(
            "atlas_feat_001", 700.456,
            [150.05, 250.10, 350.15], [1.0, 0.7, 0.3],
        ),
    ]
    atlas_mgf = tmp_path / "atlas.mgf"
    save_as_mgf(atlas, str(atlas_mgf))

    # Soil: 1 spectrum matches atlas_feat_000 (same precursor, same peaks),
    # 1 doesn't match anything in the atlas.
    soil = [
        _make_spectrum(
            "soil_01", 500.123,
            [100.05, 200.10, 300.15, 400.20], [1.0, 0.8, 0.6, 0.4],
        ),
        _make_spectrum("soil_02", 999.999, [999.0], [1.0]),
    ]
    soil_mgf = tmp_path / "soil.mgf"
    save_as_mgf(soil, str(soil_mgf))

    # Minimal SIMPER atlas: atlas_feat_000 -> Bacteria/Bacillota,
    # atlas_feat_001 -> Fungi/Ascomycota
    simper = pd.DataFrame([
        {
            "feature_id": "atlas_feat_000", "phylum": "Bacillota",
            "kingdom": "Bacteria", "direction": "enriched",
            "fold_change": 5.0, "simper_rank": 0,
        },
        {
            "feature_id": "atlas_feat_001", "phylum": "Ascomycota",
            "kingdom": "Fungi", "direction": "enriched",
            "fold_change": 3.5, "simper_rank": 0,
        },
    ])
    simper_p = tmp_path / "simper.parquet"
    simper.to_parquet(simper_p)

    out = tmp_path / "matches.parquet"
    subprocess.run([
        sys.executable, str(BIN),
        "--soil-mgf", str(soil_mgf),
        "--atlas-mgf", str(atlas_mgf),
        "--simper-atlas", str(simper_p),
        "--min-cos", "0.5",          # loose for test
        "--min-matched-peaks", "2",  # loose for test
        "--output-matches", str(out),
    ], check=True)

    df = pd.read_parquet(out)

    # Expected: exactly one row — soil_01 x atlas_feat_000 — after filters.
    assert len(df) == 1, f"expected 1 match row, got {len(df)}: {df}"
    assert df["soil_feature_id"].iloc[0] == "soil_01"
    assert df["atlas_feature_id"].iloc[0] == "atlas_feat_000"
    assert df["phylum"].iloc[0] == "Bacillota"
    assert df["kingdom"].iloc[0] == "Bacteria"
    assert df["direction"].iloc[0] == "enriched"
    assert df["cosine"].iloc[0] >= 0.5
    assert df["n_matched_peaks"].iloc[0] >= 2
    # Same precursor -> ppm diff should be tiny
    assert df["precursor_ppm_diff"].iloc[0] < 1.0

    # Full 10-column schema contract.
    expected_cols = {
        "soil_feature_id", "atlas_feature_id", "phylum", "kingdom",
        "direction", "fold_change", "simper_rank", "cosine",
        "n_matched_peaks", "precursor_ppm_diff",
    }
    assert expected_cols.issubset(set(df.columns)), (
        f"missing columns: {expected_cols - set(df.columns)}"
    )


def test_simper_match_empty_result(tmp_path):
    """When no spectra match, output should be an empty parquet with correct
    schema (no crash, no KeyError, no missing columns)."""
    soil = [_make_spectrum("soil_01", 100.0, [50.0], [1.0])]
    atlas = [_make_spectrum("atlas_01", 999.0, [50.0], [1.0])]
    soil_mgf = tmp_path / "soil.mgf"
    atlas_mgf = tmp_path / "atlas.mgf"
    save_as_mgf(soil, str(soil_mgf))
    save_as_mgf(atlas, str(atlas_mgf))

    simper = pd.DataFrame([
        {
            "feature_id": "atlas_01", "phylum": "Bacillota",
            "kingdom": "Bacteria", "direction": "enriched",
            "fold_change": 5.0, "simper_rank": 0,
        },
    ])
    simper_p = tmp_path / "simper.parquet"
    simper.to_parquet(simper_p)

    out = tmp_path / "matches.parquet"
    subprocess.run([
        sys.executable, str(BIN),
        "--soil-mgf", str(soil_mgf),
        "--atlas-mgf", str(atlas_mgf),
        "--simper-atlas", str(simper_p),
        "--output-matches", str(out),
    ], check=True)

    df = pd.read_parquet(out)
    assert len(df) == 0
    # Schema must still be present on empty output — downstream concat logic
    # depends on stable columns.
    expected_cols = {
        "soil_feature_id", "atlas_feature_id", "phylum", "kingdom",
        "direction", "fold_change", "simper_rank", "cosine",
        "n_matched_peaks", "precursor_ppm_diff",
    }
    assert expected_cols.issubset(set(df.columns)), (
        f"empty output missing columns: {expected_cols - set(df.columns)}"
    )


def test_simper_match_help_does_not_crash():
    result = subprocess.run(
        [sys.executable, str(BIN), "--help"],
        capture_output=True, text=True, check=True,
    )
    assert "--soil-mgf" in result.stdout
    assert "--atlas-mgf" in result.stdout
    assert "--simper-atlas" in result.stdout
    assert "--output-matches" in result.stdout
