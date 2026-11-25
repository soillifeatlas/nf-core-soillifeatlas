import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import application as ap  # noqa: E402


def test_anova_finds_a_pure_drought_effect():
    d = np.array([0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1], float)
    c = np.array([0, 0, 0, 1, 1, 1, 0, 0, 0, 1, 1, 1], float)
    y = 2 + 1.0 * d + np.array([.01, -.02, .01, .02, -.01, 0, .01, 0, -.01, .02, -.02, .01])
    r = ap.anova2(y, d, c)
    assert r["p_drought"] < 1e-6 and r["p_climate"] > 0.05 and r["p_interaction"] > 0.05


def test_bh_matches_hand_calculation():
    q = ap.bh(pd.Series([0.01, 0.04, 0.03, np.nan]))
    assert np.allclose(q[:3], [0.03, 0.04, 0.04]) and np.isnan(q[3])
