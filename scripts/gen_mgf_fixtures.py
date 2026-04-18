"""One-shot: generate tiny MGF + decomposition fixtures for nf-test runs.

Writes into tests/fixtures/tiny/:
  tiny_atlas.mgf                — 2 atlas spectra (one matches soil)
  tiny_soil.mgf                 — 2 soil spectra (one matches atlas_feat_000,
                                  one has no atlas counterpart)
  tiny_simper_atlas_mgf.parquet — 2-row SIMPER atlas keyed to the MGF
                                  feature_ids above. (tiny_simper_atlas.parquet
                                  in this directory is keyed to feat_000..
                                  feat_009 and is reused for the decomposition
                                  nf-tests, but it can't be reused for the
                                  SIMPER_MATCH nf-test.)
  tiny_atlas_intensity.parquet  — 10 features x 6 samples atlas intensity
                                  matrix for DECOMPOSE nf-tests.
  tiny_sample_phylum_map.csv    — maps the 6 atlas samples to 3 phyla, 2 each.
  tiny_sample_metadata.tsv      — sample treatment assignments for
                                  EVAL_TREATMENT nf-tests.

Keep the spectra intentionally small (3-4 peaks each) so the tests stay fast.
Commit the resulting files. Re-run only if the schema needs to change.
"""
from pathlib import Path

import numpy as np
import pandas as pd
from matchms import Spectrum
from matchms.exporting import save_as_mgf


FIXTURES = Path(__file__).resolve().parents[1] / "tests" / "fixtures" / "tiny"
FIXTURES.mkdir(parents=True, exist_ok=True)


def _make_spectrum(feature_id: str, precursor_mz: float,
                   peak_mz_list: list[float], peak_int_list: list[float]) -> Spectrum:
    return Spectrum(
        mz=np.array(peak_mz_list, dtype=float),
        intensities=np.array(peak_int_list, dtype=float),
        metadata={
            "precursor_mz": precursor_mz,
            "feature_id": feature_id,
            "scans": feature_id,
        },
    )


# Atlas: 2 spectra at known precursors. atlas_feat_000 is the "real" match.
# matchms.exporting.save_as_mgf appends if the file exists, which makes
# re-running this script non-idempotent. Unlink first so the output is always
# exactly what we emit here, nothing more.
(FIXTURES / "tiny_atlas.mgf").unlink(missing_ok=True)
(FIXTURES / "tiny_soil.mgf").unlink(missing_ok=True)

atlas_spectra = [
    _make_spectrum(
        "atlas_feat_000", 500.123,
        [100.05, 200.10, 300.15, 400.20], [1.0, 0.8, 0.6, 0.4],
    ),
    _make_spectrum(
        "atlas_feat_001", 700.456,
        [150.05, 250.10, 350.15], [1.0, 0.7, 0.3],
    ),
]
save_as_mgf(atlas_spectra, str(FIXTURES / "tiny_atlas.mgf"))

# Soil: 1 matches atlas_feat_000 (same precursor + peaks), 1 has no match.
soil_spectra = [
    _make_spectrum(
        "soil_01", 500.123,
        [100.05, 200.10, 300.15, 400.20], [1.0, 0.8, 0.6, 0.4],
    ),
    _make_spectrum("soil_02", 999.999, [999.0], [1.0]),
]
save_as_mgf(soil_spectra, str(FIXTURES / "tiny_soil.mgf"))

# SIMPER atlas aligned to the MGF feature_ids above — needed because
# tiny_simper_atlas.parquet uses feat_NNN keys for the pytest wrappers.
simper_for_mgf = pd.DataFrame([
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
simper_for_mgf.to_parquet(FIXTURES / "tiny_simper_atlas_mgf.parquet")

# ---------- atlas intensity + sample->phylum map ----------
# Mirrors tests/unit/test_run_decomposition.py._build_tiny_atlas_and_soil so
# the DECOMPOSE nf-test and the Python wrapper test exercise the same shape.
# Phyla MUST match the ones baked into tiny_simper_atlas.parquet, otherwise
# build_phylum_reference_array finds no matching phyla and exits.
rng = np.random.default_rng(7)
feat_ids = [f"feat_{i:03d}" for i in range(10)]

atlas_int = pd.DataFrame(
    rng.random((10, 6)) * 1e5 + 1000,
    index=feat_ids,
    columns=["A1", "A2", "B1", "B2", "C1", "C2"],
)
atlas_int.index.name = "feature_id"
atlas_int.to_parquet(FIXTURES / "tiny_atlas_intensity.parquet")

sample_phylum = pd.DataFrame({
    "sample_id": ["A1", "A2", "B1", "B2", "C1", "C2"],
    "phylum":    ["Actinomycetota"] * 2 + ["Ascomycota"] * 2 + ["Euryarchaeota"] * 2,
})
sample_phylum.to_csv(FIXTURES / "tiny_sample_phylum_map.csv", index=False)

# ---------- sample metadata for EVAL_TREATMENT ----------
# The tiny_intensity fixture has 4 soil samples S01..S04. Split them 2/2
# across a made-up `drought` contrast so the Mann-Whitney U test has
# something to compute.
sample_meta = pd.DataFrame({
    "sample_id": ["S01", "S02", "S03", "S04"],
    "drought":   ["Yes", "Yes", "No", "No"],
})
sample_meta.to_csv(FIXTURES / "tiny_sample_metadata.tsv", sep="\t", index=False)

# ---------- tiny kingdom-composition parquet for EVAL_PLAUSIBILITY ----------
# Long format, matches what run_decomposition emits: one row per
# (sample_id, kingdom), proportion_pct sums to 100 per sample.
# Use the same kingdoms present in tiny_expected_kingdom_composition.csv.
kingdom_comp = pd.DataFrame([
    {"sample_id": "SOIL_01", "kingdom": "Bacteria", "proportion_pct": 40.0},
    {"sample_id": "SOIL_01", "kingdom": "Fungi",    "proportion_pct": 25.0},
    {"sample_id": "SOIL_01", "kingdom": "Archaea",  "proportion_pct": 5.0},
    {"sample_id": "SOIL_01", "kingdom": "Plantae",  "proportion_pct": 28.0},
    {"sample_id": "SOIL_01", "kingdom": "Animalia", "proportion_pct": 1.5},
    {"sample_id": "SOIL_01", "kingdom": "Protozoa", "proportion_pct": 0.5},
    {"sample_id": "SOIL_02", "kingdom": "Bacteria", "proportion_pct": 45.0},
    {"sample_id": "SOIL_02", "kingdom": "Fungi",    "proportion_pct": 22.0},
    {"sample_id": "SOIL_02", "kingdom": "Archaea",  "proportion_pct": 4.0},
    {"sample_id": "SOIL_02", "kingdom": "Plantae",  "proportion_pct": 26.0},
    {"sample_id": "SOIL_02", "kingdom": "Animalia", "proportion_pct": 2.5},
    {"sample_id": "SOIL_02", "kingdom": "Protozoa", "proportion_pct": 0.5},
])
kingdom_comp.to_parquet(FIXTURES / "tiny_kingdom_composition.parquet")

# ---------- tiny phylum-composition parquet for EVAL_TREATMENT ----------
# One row per (sample_id, phylum). Uses the 4 sample metadata sample_ids so
# the MWU contrast has actual Yes/No splits to work with.
phylum_comp = pd.DataFrame([
    {"sample_id": "S01", "phylum": "Actinomycetota", "proportion_pct": 40.0},
    {"sample_id": "S01", "phylum": "Ascomycota",     "proportion_pct": 30.0},
    {"sample_id": "S01", "phylum": "Euryarchaeota",  "proportion_pct": 30.0},
    {"sample_id": "S02", "phylum": "Actinomycetota", "proportion_pct": 38.0},
    {"sample_id": "S02", "phylum": "Ascomycota",     "proportion_pct": 32.0},
    {"sample_id": "S02", "phylum": "Euryarchaeota",  "proportion_pct": 30.0},
    {"sample_id": "S03", "phylum": "Actinomycetota", "proportion_pct": 30.0},
    {"sample_id": "S03", "phylum": "Ascomycota",     "proportion_pct": 40.0},
    {"sample_id": "S03", "phylum": "Euryarchaeota",  "proportion_pct": 30.0},
    {"sample_id": "S04", "phylum": "Actinomycetota", "proportion_pct": 28.0},
    {"sample_id": "S04", "phylum": "Ascomycota",     "proportion_pct": 42.0},
    {"sample_id": "S04", "phylum": "Euryarchaeota",  "proportion_pct": 30.0},
])
phylum_comp.to_parquet(FIXTURES / "tiny_phylum_composition.parquet")

print(f"MGF fixtures written to: {FIXTURES}")
