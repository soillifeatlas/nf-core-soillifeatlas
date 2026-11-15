"""Unit tests for the stage 3 quality checks (soillifeatlas.process_qc)."""
from __future__ import annotations

from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from soillifeatlas.process_qc import (
    Thresholds,
    file_id,
    match_internal_standards,
    normalise_role,
    per_file_matrix,
    read_mzml_start_time,
    read_samplesheet,
    run_qc,
)

PC_MZ = 753.61338  # PC 15:0-18:1(d7) [M+H]+
FILES = ["S1", "S2", "QC1", "QC2", "QC3", "BLK"]


def _targets() -> pd.DataFrame:
    return pd.DataFrame(
        [
            {"compound": "PC_15_0_18_1_d7", "lipid_class": "PC", "polarity": "pos", "adduct": "[M+H]+",
             "mz": PC_MZ, "priority": "primary", "spike_ng": 10},
            {"compound": "PC_15_0_18_1_d7", "lipid_class": "PC", "polarity": "pos", "adduct": "[M+Na]+",
             "mz": 775.59532, "priority": "secondary", "spike_ng": 10},
            {"compound": "LPE_18_1_d7", "lipid_class": "LPE", "polarity": "pos", "adduct": "[M+H]+",
             "mz": 487.35238, "priority": "primary", "spike_ng": 10},
        ]
    )


def _full(pc_areas, pc_rts, pc_mz_shift_ppm=1.0) -> pd.DataFrame:
    """Three features: the PC standard, its sodium adduct (smaller), and a lipid."""
    rows = []
    features = [
        (1, PC_MZ * (1 + pc_mz_shift_ppm * 1e-6), 15.6, pc_areas, pc_rts),
        (2, 775.59532, 15.6, [a * 0.1 for a in pc_areas], pc_rts),
        (3, 700.5, 12.0, [1000, 1100, 1000, 1010, 990, 0], [12.0] * 6),
    ]
    for fid, mz, rt, areas, rts in features:
        r = {"id": fid, "mz": mz, "rt": rt}
        for f, a, t in zip(FILES, areas, rts):
            r[f"datafile:{f}.mzML:area"] = a
            r[f"datafile:{f}.mzML:mz"] = mz if a else np.nan
            r[f"datafile:{f}.mzML:rt"] = t if a else np.nan
        rows.append(r)
    return pd.DataFrame(rows)


def _sheet(tmp_path: Path) -> pd.DataFrame:
    p = tmp_path / "sheet.csv"
    pd.DataFrame(
        {
            "file": [f"/data/{f}.raw" for f in FILES],
            "polarity": ["pos"] * 6,
            "role": ["sample", "sample", "QC", "pooled_QC", "pooled_qc", "solvent_blank"],
        }
    ).to_csv(p, index=False)
    return read_samplesheet(p)


def test_file_id_strips_path_and_extension():
    assert file_id("/a/b/OE12-3-19POS.raw") == "OE12-3-19POS"
    assert file_id("OE12-3-19POS.mzML") == "OE12-3-19POS"


def test_role_aliases():
    assert normalise_role("QC") == "pooled_qc"
    assert normalise_role("SRM") == "nist"
    assert normalise_role("solvent blank") == "blank"


def test_samplesheet_rejects_unknown_role(tmp_path):
    p = tmp_path / "bad.csv"
    pd.DataFrame({"file": ["a.raw"], "polarity": ["pos"], "role": ["mystery"]}).to_csv(p, index=False)
    with pytest.raises(ValueError, match="Unknown role"):
        read_samplesheet(p)


def test_samplesheet_rejects_bad_polarity(tmp_path):
    p = tmp_path / "bad.csv"
    pd.DataFrame({"file": ["a.raw"], "polarity": ["positive"], "role": ["sample"]}).to_csv(p, index=False)
    with pytest.raises(ValueError, match="polarity"):
        read_samplesheet(p)


def test_per_file_matrix_uses_file_ids():
    full = _full([100] * 6, [15.6] * 6)
    area = per_file_matrix(full, "area")
    assert list(area.columns) == FILES
    assert area.loc[3, "BLK"] == 0


def test_primary_adduct_wins_over_larger_secondary():
    full = _full([100] * 6, [15.6] * 6)
    area = per_file_matrix(full, "area")
    m = match_internal_standards(full, _targets(), "pos", area).set_index("compound")
    assert m.loc["PC_15_0_18_1_d7", "adduct"] == "[M+H]+"
    assert m.loc["PC_15_0_18_1_d7", "feature_id"] == 1
    assert not m.loc["LPE_18_1_d7", "matched"]


def test_good_run_passes_is_checks(tmp_path):
    full = _full([1000, 1050, 980, 1010, 1000, 0], [15.60, 15.62, 15.61, 15.59, 15.60, 15.6])
    res = run_qc(full, _sheet(tmp_path), _targets(), "pos")
    pc = res.internal_standards.set_index("compound").loc["PC_15_0_18_1_d7"]
    assert pc["area_rsd_pct"] < 5
    assert pc["max_rt_deviation_min"] == pytest.approx(0.02, abs=1e-6)
    status = res.checks_frame().set_index("check")["status"]
    assert status["Internal standard retention time"] == "pass"
    assert status["Internal standard peak area RSD (samples and pooled QC)"] == "pass"
    assert status["Internal standard mass accuracy"] == "pass"
    # LPE target is absent, so not all standards are found
    assert status["Internal standards found"] == "warn"


def test_drifting_standard_fails(tmp_path):
    full = _full([1000, 400, 1500, 900, 1000, 0], [15.6, 15.9, 15.6, 15.6, 15.6, 15.6], pc_mz_shift_ppm=8)
    res = run_qc(full, _sheet(tmp_path), _targets(), "pos")
    status = res.checks_frame().set_index("check")["status"]
    assert status["Internal standard retention time"] == "fail"
    assert status["Internal standard peak area RSD (samples and pooled QC)"] == "fail"
    assert status["Internal standard mass accuracy"] == "fail"


def test_blank_is_excluded_from_is_rsd(tmp_path):
    # A blank with a huge IS area must not inflate the RSD.
    full = _full([1000, 1000, 1000, 1000, 1000, 50000], [15.6] * 6)
    res = run_qc(full, _sheet(tmp_path), _targets(), "pos")
    pc = res.internal_standards.set_index("compound").loc["PC_15_0_18_1_d7"]
    assert pc["area_rsd_pct"] == pytest.approx(0.0, abs=1e-9)


def test_pooled_qc_rsd_and_blank_ratio(tmp_path):
    full = _full([1000] * 5 + [900], [15.6] * 6)
    res = run_qc(full, _sheet(tmp_path), _targets(), "pos")
    f = res.features
    assert f.loc[3, "qc_rsd_pct"] < 2
    # feature 1: blank 900 vs sample median 1000
    assert f.loc[1, "blank_ratio"] == pytest.approx(0.9)


def test_injection_order_from_start_times(tmp_path):
    full = _full([100] * 6, [15.6] * 6)
    times = {f: datetime(2023, 3, 30, 10 - i, tzinfo=timezone.utc) for i, f in enumerate(FILES)}
    res = run_qc(full, _sheet(tmp_path), _targets(), "pos", start_times=times)
    assert list(res.injections["file_id"]) == list(reversed(FILES))
    assert res.injections["order_source"].iloc[0] == "mzML start time"


def test_read_mzml_start_time(tmp_path):
    p = tmp_path / "x.mzML"
    p.write_text('<mzML><run id="x" startTimeStamp="2023-03-30T00:21:47Z"></run></mzML>')
    assert read_mzml_start_time(p) == datetime(2023, 3, 30, 0, 21, 47, tzinfo=timezone.utc)
    q = tmp_path / "y.mzML"
    q.write_text("<mzML></mzML>")
    assert read_mzml_start_time(q) is None


def test_thresholds_are_configurable(tmp_path):
    full = _full([1000, 1300, 1000, 1000, 1000, 0], [15.6] * 6)
    strict = run_qc(full, _sheet(tmp_path), _targets(), "pos", thresholds=Thresholds(is_area_rsd_pct=5))
    loose = run_qc(full, _sheet(tmp_path), _targets(), "pos", thresholds=Thresholds(is_area_rsd_pct=50))
    key = "Internal standard peak area RSD (samples and pooled QC)"
    assert strict.checks_frame().set_index("check").loc[key, "status"] == "fail"
    assert loose.checks_frame().set_index("check").loc[key, "status"] == "pass"


def test_standard_injection_beats_isobaric_soil_lipid():
    """A soil lipid at the [M+H]+ mass must not be taken for the standard.

    Feature 10: soil lipid at the PC-d7 [M+H]+ m/z, large in soils, absent
    from the standard-only injection. Feature 11: the real standard, seen as
    [M+Na]+ and dominant in the standard-only injection.
    """
    files = ["S1", "S2", "ISTD"]
    rows = []
    for fid, mz, rt, areas in [
        (10, PC_MZ, 9.1, [5e5, 6e5, 0]),
        (11, 775.59532, 15.6, [1e4, 1.1e4, 8e5]),
    ]:
        r = {"id": fid, "mz": mz, "rt": rt}
        for f, a in zip(files, areas):
            r[f"datafile:{f}.mzML:area"] = a
            r[f"datafile:{f}.mzML:mz"] = mz
            r[f"datafile:{f}.mzML:rt"] = rt
        rows.append(r)
    full = pd.DataFrame(rows)
    area = per_file_matrix(full, "area")
    without = match_internal_standards(full, _targets(), "pos", area).set_index("compound")
    assert without.loc["PC_15_0_18_1_d7", "feature_id"] == 10  # the trap, without a reference
    with_ref = match_internal_standards(full, _targets(), "pos", area, reference_files=["ISTD"]).set_index("compound")
    assert with_ref.loc["PC_15_0_18_1_d7", "feature_id"] == 11
    assert with_ref.loc["PC_15_0_18_1_d7", "adduct"] == "[M+Na]+"
    assert with_ref.loc["PC_15_0_18_1_d7", "selected_by"] == "standard injection"


def test_custom_standards_table_with_minimal_columns(tmp_path):
    from soillifeatlas.process_qc import read_targets

    p = tmp_path / "my_standards.csv"
    pd.DataFrame({"compound": ["SPLASH_PC"], "polarity": ["POS"], "adduct": ["[M+H]+"], "mz": [PC_MZ]}).to_csv(p, index=False)
    tg = read_targets(p)
    assert tg.loc[0, "priority"] == "primary"
    assert tg.loc[0, "lipid_class"] == "SPLASH_PC"
    assert tg.loc[0, "polarity"] == "pos"


def test_subset_of_shipped_equisplash_table():
    from soillifeatlas.process_qc import read_targets

    shipped = Path(__file__).resolve().parents[2] / "assets" / "equisplash_targets.csv"
    tg = read_targets(shipped, ["PC_15_0_18_1_d7", "LPE_18_1_d7"])
    assert set(tg["compound"]) == {"PC_15_0_18_1_d7", "LPE_18_1_d7"}
    with pytest.raises(ValueError, match="not in the table"):
        read_targets(shipped, ["PC_99_9"])


def test_expected_rt_selects_the_right_feature():
    files = ["S1", "S2"]
    rows = []
    for fid, rt, areas in [(20, 9.1, [9e5, 9e5]), (21, 15.6, [1e5, 1e5])]:
        r = {"id": fid, "mz": PC_MZ, "rt": rt}
        for f, a in zip(files, areas):
            r[f"datafile:{f}.mzML:area"] = a
            r[f"datafile:{f}.mzML:mz"] = PC_MZ
            r[f"datafile:{f}.mzML:rt"] = rt
        rows.append(r)
    full = pd.DataFrame(rows)
    area = per_file_matrix(full, "area")
    tg = _targets()
    tg["rt"] = [15.55, 15.55, np.nan]
    m = match_internal_standards(full, tg, "pos", area, rt_window=0.3).set_index("compound")
    assert m.loc["PC_15_0_18_1_d7", "feature_id"] == 21
    assert m.loc["PC_15_0_18_1_d7", "selected_by"] == "expected RT"
