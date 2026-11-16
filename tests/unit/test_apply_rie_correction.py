"""Unit test for bin/apply_rie_correction.py.

Exercises the CLI via subprocess and verifies:
  * output parquet is produced with the expected shape
  * the RIE floor caps amplification for very-low-RIE classes (MG, RIE=0.0005)
    so we do NOT get a 2000x blowup
"""
import subprocess
import sys
from pathlib import Path

import pandas as pd

FIXTURES = Path(__file__).resolve().parents[1] / "fixtures" / "tiny"
REPO_ROOT = Path(__file__).resolve().parents[2]
BIN = REPO_ROOT / "bin" / "apply_rie_correction.py"


def test_apply_rie_correction_with_floor(tmp_path):
    out_file = tmp_path / "rie_corrected.parquet"
    subprocess.run(
        [
            sys.executable, str(BIN),
            "--intensity", str(FIXTURES / "tiny_intensity.parquet"),
            "--rie-table", str(FIXTURES / "tiny_rie_table.csv"),
            "--annotation", str(FIXTURES / "tiny_annotation.csv"),
            "--rie-floor", "0.20",
            "--output", str(out_file),
        ],
        check=True,
    )
    assert out_file.exists()

    raw = pd.read_parquet(FIXTURES / "tiny_intensity.parquet")
    corrected = pd.read_parquet(out_file)

    assert corrected.shape == (10, 4)
    assert list(corrected.index) == list(raw.index)
    assert list(corrected.columns) == list(raw.columns)

    # feat_002 is class MG with RIE=0.0005. Without floor, correction would be
    # 1/0.0005 = 2000x amplification. With floor=0.20, max amplification is
    # 1/0.20 = 5x. So corrected_MG_row / raw_MG_row should be ~5 (exactly 5.0
    # given corrected = intensity / clip(RIE, floor, ceiling) = intensity / 0.20).
    raw_mg = raw.loc["feat_002"]
    corrected_mg = corrected.loc["feat_002"]
    ratio = (corrected_mg / raw_mg).values
    # Exact: 1/0.20 = 5.0
    assert all(abs(r - 5.0) < 1e-9 for r in ratio), (
        f"MG row (RIE=0.0005) should be amplified by exactly 1/floor = 5.0 "
        f"with floor=0.20, got ratios {ratio}"
    )

    # feat_001 (class PC, RIE=1.0) should be essentially unchanged (divided by 1.0)
    raw_pc = raw.loc["feat_001"]
    corrected_pc = corrected.loc["feat_001"]
    ratio_pc = (corrected_pc / raw_pc).values
    assert all(abs(r - 1.0) < 1e-9 for r in ratio_pc), (
        f"PC row (RIE=1.0) should be unchanged, got ratios {ratio_pc}"
    )


def test_apply_rie_correction_help_does_not_crash():
    result = subprocess.run(
        [sys.executable, str(BIN), "--help"],
        capture_output=True,
        text=True,
        check=True,
    )
    assert "--intensity" in result.stdout
    assert "--rie-table" in result.stdout
    assert "--annotation" in result.stdout
    assert "--rie-floor" in result.stdout
