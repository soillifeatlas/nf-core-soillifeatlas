"""Unit test for bin/apply_ref_quality_filter.py.

Verifies Layer 5 (reference-quality filter): Archaea SIMPER atlas rows are
retained only when their feature_id appears in the ArchLips-validated set.
Non-Archaea rows are untouched.
"""
import subprocess
import sys
from pathlib import Path

import pandas as pd

FIXTURES = Path(__file__).resolve().parents[1] / "fixtures" / "tiny"
REPO_ROOT = Path(__file__).resolve().parents[2]
BIN = REPO_ROOT / "bin" / "apply_ref_quality_filter.py"


def test_ref_quality_filter_keeps_only_archlips_archaea(tmp_path):
    out_file = tmp_path / "filtered_atlas.parquet"
    subprocess.run(
        [
            sys.executable, str(BIN),
            "--simper-atlas", str(FIXTURES / "tiny_simper_atlas.parquet"),
            "--archlips-validated", str(FIXTURES / "tiny_archlips_validated.csv"),
            "--output", str(out_file),
        ],
        check=True,
    )
    assert out_file.exists()

    raw = pd.read_parquet(FIXTURES / "tiny_simper_atlas.parquet")
    filtered = pd.read_parquet(out_file)

    validated = set(pd.read_csv(FIXTURES / "tiny_archlips_validated.csv")["feature_id"])
    # Non-Archaea rows are unchanged
    raw_non_archaea = raw[raw["kingdom"] != "Archaea"].reset_index(drop=True)
    filt_non_archaea = filtered[filtered["kingdom"] != "Archaea"].reset_index(drop=True)
    pd.testing.assert_frame_equal(raw_non_archaea, filt_non_archaea)

    # Every Archaea row in the output has a feature_id in the validated set
    filt_archaea = filtered[filtered["kingdom"] == "Archaea"]
    assert not filt_archaea.empty, "expected at least one validated Archaea row to survive"
    assert set(filt_archaea["feature_id"]) <= validated, (
        f"filtered Archaea contains non-validated features: "
        f"{set(filt_archaea['feature_id']) - validated}"
    )

    # All Archaea rows that WERE in the validated set must be kept
    raw_valid_archaea = raw[(raw["kingdom"] == "Archaea") & raw["feature_id"].isin(validated)]
    assert len(filt_archaea) == len(raw_valid_archaea), (
        f"expected {len(raw_valid_archaea)} validated Archaea rows kept, "
        f"got {len(filt_archaea)}"
    )

    # Sanity: non-validated Archaea rows were actually dropped (the fixture has
    # 10 Archaea rows and 3 validated feature_ids => 7 should be dropped)
    assert len(raw[raw["kingdom"] == "Archaea"]) == 10
    assert len(filt_archaea) == 3


def test_ref_quality_filter_help_does_not_crash():
    result = subprocess.run(
        [sys.executable, str(BIN), "--help"],
        capture_output=True,
        text=True,
        check=True,
    )
    assert "--simper-atlas" in result.stdout
    assert "--archlips-validated" in result.stdout
