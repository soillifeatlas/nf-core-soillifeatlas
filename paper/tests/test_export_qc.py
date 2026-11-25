import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import export_qc as ex  # noqa: E402


def test_gapfilled_values_enter_only_the_second_rsd():
    area = pd.DataFrame({"q1": [100.0, 100.0], "q2": [100.0, 100.0], "q3": [100.0, 400.0]}, index=[1, 2])
    det = pd.DataFrame({"q1": [True, True], "q2": [True, True], "q3": [True, False]}, index=[1, 2])
    per, summ = ex.rsd_tables(area, det, ["q1", "q2", "q3"], 0.6, 30)
    assert per.loc[1, "rsd_peak_picked_pct"] == 0
    assert per.loc[2, "rsd_peak_picked_pct"] == 0          # q3 gap-filled, left out
    assert per.loc[2, "rsd_with_gapfill_pct"] > 30          # q3 counted
    assert summ.set_index("mode").loc["peak_picked_only", "n_eligible"] == 2


def test_is_amount_against_calibration():
    rows = [{"compound": "PC", "file_id": f"OE12-1ESTD-A{i}POS", "role": "standard_cal", "detected": True,
             "area": 1000 * ex.ESTD1_NG_PER_UL[f"A{i}"]} for i in range(1, 8)]
    rows.append({"compound": "PC", "file_id": "s1", "role": "sample", "detected": True, "area": 20.0})
    t = ex.is_vs_calibration(pd.DataFrame(rows), recon_ul=100, spike_ng=10).iloc[0]
    assert np.isclose(t["soil_median_ng_per_ul"], 0.02)
    assert t["status"] == "below"
    assert np.isclose(t["R_app_median_pct"], 20.0)
