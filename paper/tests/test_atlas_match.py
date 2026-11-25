import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import atlas_match as am  # noqa: E402


def test_candidates_respect_ppm_window():
    lib = np.array([500.0, 500.004, 500.006, 600.0])
    assert list(am.candidates(500.0, lib, 10)) == [0, 1]   # 10 ppm of 500 = 0.005


def test_detection_rate_uses_units_per_phylum():
    pres = pd.DataFrame({"feature_id": ["a", "a", "b", "b", "b"],
                         "unit": ["u1", "u2", "u1", "u3", "u4"],
                         "phylum": ["P1", "P1", "P1", "P2", "P2"]})
    d = am.detection_by_phylum(pres, ["a"]).set_index("phylum")
    assert d.loc["P1", "detection_rate"] == 1.0          # a in u1, u2 of P1's two units
    assert "P2" not in d.index
