"""Unit test for bin/eval_treatment_effects.py.

TDD: written before the wrapper. Builds a synthetic phylum-level composition
parquet with designed-in effects (Actinomycetota enriched under drought),
then verifies the wrapper's Mann-Whitney U output reflects that direction.

Sample design: 4 treatments x 3 samples each = 12 samples total, split across
  - Ambient_Control (drought=No, warming=No)
  - Ambient_Drought (drought=Yes, warming=No)
  - Future_Control  (drought=No, warming=Yes)
  - Future_Drought  (drought=Yes, warming=Yes)

Three phyla: Actinomycetota, Ascomycota, Euryarchaeota.
Contracts asserted:
  * `Actinomycetota x drought` row has direction == "+"
  * All 6 (phylum x contrast) rows present
  * p_value column is numeric
"""
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd


REPO_ROOT = Path(__file__).resolve().parents[2]
BIN = REPO_ROOT / "bin" / "eval_treatment_effects.py"


def _build_composition_and_metadata(tmp_path: Path) -> tuple[Path, Path]:
    """12 soil samples, 3 phyla, designed with drought-enriched Actinomycetota."""
    rng = np.random.default_rng(11)

    treatments = [
        # (sample_id, treatment, drought, warming)
        ("SOIL_01", "Ambient_Control", "No",  "No"),
        ("SOIL_02", "Ambient_Control", "No",  "No"),
        ("SOIL_03", "Ambient_Control", "No",  "No"),
        ("SOIL_04", "Ambient_Drought", "Yes", "No"),
        ("SOIL_05", "Ambient_Drought", "Yes", "No"),
        ("SOIL_06", "Ambient_Drought", "Yes", "No"),
        ("SOIL_07", "Future_Control",  "No",  "Yes"),
        ("SOIL_08", "Future_Control",  "No",  "Yes"),
        ("SOIL_09", "Future_Control",  "No",  "Yes"),
        ("SOIL_10", "Future_Drought",  "Yes", "Yes"),
        ("SOIL_11", "Future_Drought",  "Yes", "Yes"),
        ("SOIL_12", "Future_Drought",  "Yes", "Yes"),
    ]
    meta = pd.DataFrame(treatments, columns=["sample_id", "treatment", "drought", "warming"])
    meta_p = tmp_path / "sample_metadata.tsv"
    meta.to_csv(meta_p, sep="\t", index=False)

    # Designed proportions: drought -> Actinomycetota UP (shift means), other
    # phyla take the complement. Small jitter via rng so the ranksum test
    # has variance but a clean direction.
    rows = []
    for sid, _treat, drought, _warm in treatments:
        acti = 60.0 + rng.normal(0, 2.0) if drought == "Yes" else 30.0 + rng.normal(0, 2.0)
        asco = 25.0 + rng.normal(0, 2.0) if drought == "Yes" else 40.0 + rng.normal(0, 2.0)
        eury = max(0.0, 100.0 - acti - asco)
        # Renormalise so each sample sums to exactly 100 despite the jitter
        total = acti + asco + eury
        acti, asco, eury = (100 * v / total for v in (acti, asco, eury))
        rows += [
            {"sample_id": sid, "phylum": "Actinomycetota", "proportion_pct": acti},
            {"sample_id": sid, "phylum": "Ascomycota",     "proportion_pct": asco},
            {"sample_id": sid, "phylum": "Euryarchaeota",  "proportion_pct": eury},
        ]
    comp = pd.DataFrame(rows)
    comp_p = tmp_path / "composition_phylum.parquet"
    comp.to_parquet(comp_p)
    return comp_p, meta_p


def test_eval_treatment_effects_preserves_actinomycetota_drought_direction(tmp_path):
    comp_p, meta_p = _build_composition_and_metadata(tmp_path)
    out = tmp_path / "treatment_effects.tsv"

    subprocess.run(
        [
            sys.executable, str(BIN),
            "--composition-phylum", str(comp_p),
            "--sample-metadata", str(meta_p),
            "--contrast", "drought", "warming",
            "--output", str(out),
        ],
        check=True,
    )
    assert out.exists(), "wrapper did not produce treatment_effects TSV"

    df = pd.read_csv(out, sep="\t")
    required_cols = {"phylum", "contrast", "direction", "p_value",
                     "median_yes", "median_no", "n_yes", "n_no"}
    assert required_cols.issubset(df.columns), (
        f"missing cols: {required_cols - set(df.columns)}"
    )

    # Every (phylum x contrast) combination present (3 phyla x 2 contrasts = 6 rows)
    assert len(df) == 6, f"expected 6 rows, got {len(df)}:\n{df}"
    pairs = set(zip(df["phylum"], df["contrast"]))
    assert pairs == {
        ("Actinomycetota", "drought"),
        ("Ascomycota",     "drought"),
        ("Euryarchaeota",  "drought"),
        ("Actinomycetota", "warming"),
        ("Ascomycota",     "warming"),
        ("Euryarchaeota",  "warming"),
    }

    # The target assertion: Actinomycetota UP under drought
    acti_drought = df[(df["phylum"] == "Actinomycetota") & (df["contrast"] == "drought")]
    assert len(acti_drought) == 1
    assert acti_drought["direction"].iloc[0] == "+", (
        f"expected Actinomycetota drought direction '+', got {acti_drought['direction'].iloc[0]!r}"
    )

    # p_value numeric
    assert pd.api.types.is_numeric_dtype(df["p_value"]), (
        f"p_value should be numeric, got dtype {df['p_value'].dtype}"
    )
    # n_yes / n_no should be 6 each per contrast (half the 12 samples)
    for _, row in df.iterrows():
        assert row["n_yes"] == 6
        assert row["n_no"] == 6


def test_eval_treatment_effects_help_does_not_crash():
    result = subprocess.run(
        [sys.executable, str(BIN), "--help"],
        capture_output=True, text=True, check=True,
    )
    assert "--composition-phylum" in result.stdout
    assert "--sample-metadata" in result.stdout
    assert "--contrast" in result.stdout
    assert "--output" in result.stdout
