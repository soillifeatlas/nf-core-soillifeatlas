"""Unit test for bin/apply_is_scaling.py.

TDD: written before the wrapper. Exercises the CLI end-to-end via subprocess
so we verify both arg parsing and the underlying corrections.apply_IS_normalization
call.
"""
import subprocess
import sys
from pathlib import Path

import pandas as pd

FIXTURES = Path(__file__).resolve().parents[1] / "fixtures" / "tiny"
REPO_ROOT = Path(__file__).resolve().parents[2]
BIN = REPO_ROOT / "bin" / "apply_is_scaling.py"


def test_apply_is_scaling_produces_scaled_parquet(tmp_path):
    out_file = tmp_path / "corrected.parquet"
    subprocess.run(
        [
            sys.executable, str(BIN),
            "--intensity", str(FIXTURES / "tiny_intensity.parquet"),
            "--is-features", str(FIXTURES / "tiny_is_features.csv"),
            "--is-spiked-pmol", "100",
            "--is-reference", "LPE_18d7",
            "--output", str(out_file),
        ],
        check=True,
    )
    assert out_file.exists(), "wrapper did not produce an output parquet"

    raw = pd.read_parquet(FIXTURES / "tiny_intensity.parquet")
    corrected = pd.read_parquet(out_file)

    # Shape preserved (10 features x 4 samples)
    assert corrected.shape == (10, 4)
    # Column + index order preserved
    assert list(corrected.columns) == list(raw.columns)
    assert list(corrected.index) == list(raw.index)
    # IS normalization multiplies by IS_spiked_pmol / IS_signal.
    # With random intensities near 1e6 and spiked_pmol=100, corrected values land
    # in roughly the 0..2000 range (not in the raw 1e5-1e6 range).
    assert (corrected.values >= 0).all(), "IS-normalized values must be non-negative"
    assert corrected.values.max() < raw.values.max(), (
        "IS normalization should shrink the dynamic range when spiked_pmol=100 "
        f"(raw_max={raw.values.max():.2e}, corrected_max={corrected.values.max():.2e})"
    )


def test_apply_is_scaling_help_does_not_crash():
    result = subprocess.run(
        [sys.executable, str(BIN), "--help"],
        capture_output=True,
        text=True,
        check=True,
    )
    assert "--intensity" in result.stdout
    assert "--is-features" in result.stdout
    assert "--is-reference" in result.stdout
