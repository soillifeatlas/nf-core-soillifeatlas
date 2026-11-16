"""Unit test for bin/diagnose_top_features.py (v0.1).

TDD: written before the wrapper. The v0.1 scope is intentionally minimal —
rank top-N features per phylum from the SIMPER atlas and emit them with their
mean intensity across samples for eyeball inspection. The full "RIE
amplification detector" from analysis-19 is deferred to v0.2.
"""
import subprocess
import sys
from pathlib import Path

import pandas as pd


FIXTURES = Path(__file__).resolve().parents[1] / "fixtures" / "tiny"
REPO_ROOT = Path(__file__).resolve().parents[2]
BIN = REPO_ROOT / "bin" / "diagnose_top_features.py"


def test_diagnose_top_features_emits_ranked_tsv(tmp_path):
    out = tmp_path / "diag.tsv"
    subprocess.run(
        [
            sys.executable, str(BIN),
            "--corrected-intensity", str(FIXTURES / "tiny_intensity.parquet"),
            "--simper-atlas", str(FIXTURES / "tiny_simper_atlas.parquet"),
            "--top-n", "5",
            "--output", str(out),
        ],
        check=True,
    )
    assert out.exists(), "wrapper did not produce diagnostic TSV"

    df = pd.read_csv(out, sep="\t")

    # All 3 fixture phyla represented
    assert set(df["phylum"].unique()) >= {"Actinomycetota", "Ascomycota", "Euryarchaeota"}

    # 3 phyla x 5 features each = 15 rows
    assert len(df) == 15, f"expected 15 rows, got {len(df)}:\n{df}"

    # Contract columns
    assert "feature_id" in df.columns
    assert "rank" in df.columns
    assert "phylum" in df.columns
    assert "mean_intensity_across_samples" in df.columns

    # Ranks are 0..top-n-1 per phylum (ascending simper_rank)
    for phylum, sub in df.groupby("phylum"):
        ranks = sorted(sub["rank"].tolist())
        assert ranks == [0, 1, 2, 3, 4], (
            f"phylum {phylum} ranks are {ranks}, expected [0,1,2,3,4]"
        )


def test_diagnose_top_features_respects_default_top_n(tmp_path):
    """Default --top-n=10; with 10 features per phylum in the fixture we
    should get 30 rows (3 phyla x 10 features)."""
    out = tmp_path / "diag_default.tsv"
    subprocess.run(
        [
            sys.executable, str(BIN),
            "--corrected-intensity", str(FIXTURES / "tiny_intensity.parquet"),
            "--simper-atlas", str(FIXTURES / "tiny_simper_atlas.parquet"),
            "--output", str(out),
        ],
        check=True,
    )
    df = pd.read_csv(out, sep="\t")
    assert len(df) == 30, f"expected 30 rows at default top-n=10, got {len(df)}"


def test_diagnose_top_features_help_does_not_crash():
    result = subprocess.run(
        [sys.executable, str(BIN), "--help"],
        capture_output=True, text=True, check=True,
    )
    assert "--corrected-intensity" in result.stdout
    assert "--simper-atlas" in result.stdout
    assert "--top-n" in result.stdout
    assert "--output" in result.stdout
