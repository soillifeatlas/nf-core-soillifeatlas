"""One-shot: generate tiny MGF fixtures for nf-test runs of SIMPER_MATCH.

Writes:
  tests/fixtures/tiny/tiny_atlas.mgf        — 2 atlas spectra (one matches soil)
  tests/fixtures/tiny/tiny_soil.mgf         — 2 soil spectra (one matches
                                              atlas_feat_000, one has no atlas
                                              counterpart)
  tests/fixtures/tiny/tiny_simper_atlas_mgf.parquet
                                            — a 2-row SIMPER atlas keyed to
                                              the MGF feature_ids above. The
                                              existing tiny_simper_atlas.parquet
                                              is keyed to feat_000..feat_009
                                              (for bin/ pytest wrappers) and
                                              cannot be reused here.

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

print(f"MGF fixtures written to: {FIXTURES}")
