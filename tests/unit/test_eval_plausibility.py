"""Unit test for bin/eval_plausibility.py.

TDD: written before the wrapper. Exercises the CLI end-to-end via subprocess
so we verify both arg parsing and the underlying
soillifeatlas.evaluation.plausibility_score call.

Contract for downstream CI: the emitted TSV must carry columns
    method, bc_vs_expected
(aliased from the upstream `bc_distance_from_midpoint` key) because the
nightly CI gates on ``plaus[plaus.method == "fc_weighted_bc"]["bc_vs_expected"]``.
"""
import subprocess
import sys
from pathlib import Path

import pandas as pd


FIXTURES = Path(__file__).resolve().parents[1] / "fixtures" / "tiny"
REPO_ROOT = Path(__file__).resolve().parents[2]
BIN = REPO_ROOT / "bin" / "eval_plausibility.py"


def _build_plausible_composition_parquet(tmp_path: Path) -> Path:
    """Single soil-sample, grassland-ish kingdom composition.

    Numbers picked to sit near the midpoint of each expected range so the
    resulting BC distance is small (well under 0.5), which the assertion
    below depends on.
    """
    comp = pd.DataFrame([
        {"sample_id": "S01", "kingdom": "Bacteria", "proportion_pct": 36.0},
        {"sample_id": "S01", "kingdom": "Fungi",    "proportion_pct": 28.0},
        {"sample_id": "S01", "kingdom": "Plantae",  "proportion_pct": 20.0},
        {"sample_id": "S01", "kingdom": "Archaea",  "proportion_pct": 10.0},
        {"sample_id": "S01", "kingdom": "Animalia", "proportion_pct": 3.0},
        {"sample_id": "S01", "kingdom": "Protozoa", "proportion_pct": 3.0},
    ])
    comp_p = tmp_path / "comp.parquet"
    comp.to_parquet(comp_p)
    return comp_p


def test_eval_plausibility_emits_tsv_with_expected_columns(tmp_path):
    comp_p = _build_plausible_composition_parquet(tmp_path)
    out = tmp_path / "plaus.tsv"

    subprocess.run(
        [
            sys.executable, str(BIN),
            "--composition-kingdom", str(comp_p),
            "--expected-ref", str(FIXTURES / "tiny_expected_kingdom_composition.csv"),
            "--method", "fc_weighted_bc",
            "--output", str(out),
        ],
        check=True,
    )
    assert out.exists(), "wrapper did not produce plausibility TSV"

    df = pd.read_csv(out, sep="\t")
    assert "method" in df.columns, f"missing 'method' column, got {list(df.columns)}"
    assert "bc_vs_expected" in df.columns, (
        f"missing 'bc_vs_expected' column (alias for bc_distance_from_midpoint), got {list(df.columns)}"
    )
    assert "in_range_fraction" in df.columns
    assert "n_kingdoms_scored" in df.columns
    assert "inflation_score" in df.columns

    # Exactly one row per method invocation
    assert len(df) == 1, f"expected 1 row, got {len(df)}"
    assert df["method"].iloc[0] == "fc_weighted_bc"

    # Our synthetic composition is close to the midpoints by construction — BC
    # should comfortably clear the loose ceiling the CI cares about (0.15).
    assert df["bc_vs_expected"].iloc[0] < 0.5, (
        f"BC vs expected should be small for near-midpoint comp, got {df['bc_vs_expected'].iloc[0]}"
    )

    # Per-kingdom deviations flattened as deviation_<kingdom> columns
    dev_cols = [c for c in df.columns if c.startswith("deviation_")]
    assert len(dev_cols) >= 1, (
        f"expected at least one deviation_<kingdom> column, got {list(df.columns)}"
    )


def test_eval_plausibility_help_does_not_crash():
    result = subprocess.run(
        [sys.executable, str(BIN), "--help"],
        capture_output=True, text=True, check=True,
    )
    assert "--composition-kingdom" in result.stdout
    assert "--expected-ref" in result.stdout
    assert "--method" in result.stdout
    assert "--output" in result.stdout
