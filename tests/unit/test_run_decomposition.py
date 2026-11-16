"""Unit test for bin/run_decomposition.py.

TDD: written before the wrapper. Dispatches into each of the four
soillifeatlas.decomposition methods (nnls, std_bc, enriched_bc, fc_weighted_bc)
via subprocess, one invocation per method — mirroring how Nextflow's `each`
operator will fan out the work in the pipeline.

Each method must emit a long-format parquet with columns
(sample_id, kingdom, proportion_pct) summing to ~100 per sample.
"""
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest


FIXTURES = Path(__file__).resolve().parents[1] / "fixtures" / "tiny"
REPO_ROOT = Path(__file__).resolve().parents[2]
BIN = REPO_ROOT / "bin" / "run_decomposition.py"


def _build_tiny_atlas_and_soil(tmp_path):
    """Build intensity atlas + sample->phylum map + soil matrix.

    Atlas: 10 features x 6 samples, samples split 2/2/2 across 3 phyla that
           match the phyla present in the committed tiny SIMPER atlas fixture.
    Soil:  10 features x 3 samples.
    """
    rng = np.random.default_rng(7)
    feat = [f"feat_{i:03d}" for i in range(10)]

    atlas = pd.DataFrame(
        rng.random((10, 6)) * 1e5 + 1000,
        index=feat,
        columns=["A1", "A2", "B1", "B2", "C1", "C2"],
    )
    atlas.index.name = "feature_id"
    atlas_parq = tmp_path / "atlas.parquet"
    atlas.to_parquet(atlas_parq)

    sample_phylum = pd.DataFrame({
        "sample_id": ["A1", "A2", "B1", "B2", "C1", "C2"],
        "phylum":    ["Actinomycetota"] * 2 + ["Ascomycota"] * 2 + ["Euryarchaeota"] * 2,
    })
    sp_csv = tmp_path / "sample_phylum.csv"
    sample_phylum.to_csv(sp_csv, index=False)

    soil = pd.DataFrame(
        rng.random((10, 3)) * 1e5 + 500,
        index=feat,
        columns=["SOIL_01", "SOIL_02", "SOIL_03"],
    )
    soil.index.name = "feature_id"
    soil_parq = tmp_path / "soil.parquet"
    soil.to_parquet(soil_parq)

    return atlas_parq, sp_csv, soil_parq


@pytest.mark.parametrize("method", ["nnls", "std_bc", "enriched_bc", "fc_weighted_bc"])
def test_run_decomposition_all_methods(tmp_path, method):
    atlas_p, sp_p, soil_p = _build_tiny_atlas_and_soil(tmp_path)
    simper = FIXTURES / "tiny_simper_atlas.parquet"
    out = tmp_path / f"comp_{method}.parquet"
    out_phylum = tmp_path / f"comp_{method}_phylum.parquet"

    subprocess.run(
        [
            sys.executable, str(BIN),
            "--soil-intensity", str(soil_p),
            "--atlas-intensity", str(atlas_p),
            "--simper-atlas", str(simper),
            "--sample-phylum-map", str(sp_p),
            "--method", method,
            "--output", str(out),
            "--output-phylum", str(out_phylum),
        ],
        check=True,
    )
    assert out.exists(), f"wrapper did not produce output parquet for method={method}"
    df = pd.read_parquet(out)

    # Contract: long-format (sample x kingdom x proportion_pct)
    assert {"sample_id", "kingdom", "proportion_pct"}.issubset(df.columns), (
        f"missing expected columns, got {list(df.columns)}"
    )

    # All three soil samples should appear
    assert set(df["sample_id"]) == {"SOIL_01", "SOIL_02", "SOIL_03"}, (
        f"expected all 3 soil samples, got {set(df['sample_id'])}"
    )

    # Proportions non-negative and each sample sums to ~100
    for s, sub in df.groupby("sample_id"):
        assert (sub["proportion_pct"] >= 0).all(), f"negative proportion in sample {s}"
        total = sub["proportion_pct"].sum()
        assert 99.0 <= total <= 101.0, f"sample {s} proportions sum to {total}, expected ~100"

    # --- Phylum-level output ------------------------------------------------
    assert out_phylum.exists(), f"wrapper did not produce phylum parquet for method={method}"
    pdf = pd.read_parquet(out_phylum)

    # Contract: long-format (sample x phylum x proportion_pct)
    assert set(pdf.columns) == {"sample_id", "phylum", "proportion_pct"}, (
        f"unexpected phylum columns, got {list(pdf.columns)}"
    )
    # All three soil samples should appear
    assert set(pdf["sample_id"]) == {"SOIL_01", "SOIL_02", "SOIL_03"}, (
        f"expected all 3 soil samples in phylum output, got {set(pdf['sample_id'])}"
    )
    # Proportions non-negative and each sample sums to ~100
    for s, sub in pdf.groupby("sample_id"):
        assert (sub["proportion_pct"] >= 0).all(), (
            f"negative proportion in sample {s} (phylum output)"
        )
        total = sub["proportion_pct"].sum()
        assert 99.0 <= total <= 101.0, (
            f"phylum sample {s} proportions sum to {total}, expected ~100"
        )


def test_run_decomposition_help_does_not_crash():
    result = subprocess.run(
        [sys.executable, str(BIN), "--help"],
        capture_output=True,
        text=True,
        check=True,
    )
    assert "--soil-intensity" in result.stdout
    assert "--atlas-intensity" in result.stdout
    assert "--method" in result.stdout


def test_run_decomposition_rejects_bad_method(tmp_path):
    atlas_p, sp_p, soil_p = _build_tiny_atlas_and_soil(tmp_path)
    simper = FIXTURES / "tiny_simper_atlas.parquet"
    out = tmp_path / "comp_bogus.parquet"

    result = subprocess.run(
        [
            sys.executable, str(BIN),
            "--soil-intensity", str(soil_p),
            "--atlas-intensity", str(atlas_p),
            "--simper-atlas", str(simper),
            "--sample-phylum-map", str(sp_p),
            "--method", "totally_made_up",
            "--output", str(out),
        ],
        capture_output=True,
        text=True,
    )
    assert result.returncode != 0, "expected non-zero exit for invalid --method"
    assert not out.exists(), "output should not be created when --method is invalid"
