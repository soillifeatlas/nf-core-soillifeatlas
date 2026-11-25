"""Unit tests for stage 5 Q1-Q4 concentration estimates (soillifeatlas.quantify)."""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from soillifeatlas.quantify import assign_routes, class_rie, estimate, range_flag


@pytest.fixture
def standards():
    # two standards: PC (class-matched for PC analytes) and TG
    return pd.DataFrame({"compound": ["PC_d7", "TG_d7"], "lipid_class": ["PC", "TG"],
                         "rt": [15.0, 25.0], "adduct": ["[M+H]+", "[M+NH4]+"]})


@pytest.fixture
def rie():
    return pd.DataFrame({"class": ["PC", "PC", "TG", "DG", "DG"],
                         "adduct": ["M+H", "M+H", "M+NH4", "M+NH4", "M+H"],
                         "RIE_LPE": [40.0, 60.0, 10.0, 5.0, 0.03]})


def test_class_rie_is_the_mean_per_class_and_adduct(rie):
    r = class_rie(rie)
    assert r[("PC", "M+H")] == 50.0
    assert r[("DG", "M+NH4")] == 5.0


def test_routes_follow_the_response_evidence(standards):
    ann = pd.DataFrame({
        "feature_id": [1, 2, 3, 4],
        "lipid_name": ["PC 34:1", "DG 34:1", "DG", "unknown"],
        "lipid_class": ["PC", "DG", "DG", None],
        "confidence": ["Gold", "Silver", "Bronze", None],
        "rt": [15.2, 21.0, 21.0, 10.0],
    })
    r = assign_routes(ann, standards, predicted=None).set_index("feature_id")
    assert r.loc[1, "quant_route"] == "Q1" and r.loc[1, "is_compound"] == "PC_d7"
    # no class-matched standard and no lipid-specific prediction: class-level response
    assert r.loc[2, "quant_route"] == "Q3"
    assert r.loc[3, "quant_route"] == "Q3"
    assert r.loc[2, "is_compound"] == "TG_d7"        # nearest standard in RT
    assert r.loc[4, "quant_route"] == "Q4"


def test_q2_needs_a_lipid_specific_prediction(standards):
    ann = pd.DataFrame({"feature_id": [2, 3], "lipid_name": ["DG 34:1", "DG"], "lipid_class": ["DG", "DG"],
                        "confidence": ["Silver", "Bronze"], "rt": [20.0, 21.0]})
    pred = pd.DataFrame({"lipid_name": ["DG 34:1", "DG"], "RIE_pred": [2.0, 9.9]})
    r = assign_routes(ann, standards, predicted=pred).set_index("feature_id")
    assert r.loc[2, "quant_route"] == "Q2"
    assert r.loc[3, "quant_route"] == "Q3"           # class-only names never get a lipid-specific value


def test_range_flag():
    assert range_flag(5.0, 10.0, 100.0) == "below"
    assert range_flag(50.0, 10.0, 100.0) == "within"
    assert range_flag(500.0, 10.0, 100.0) == "above"
    assert range_flag(np.nan, 10.0, 100.0) == "is_not_detected"


def test_estimate_applies_the_ratio_formula_without_extra_recovery(standards, rie):
    routes = pd.DataFrame({
        "feature_id": [1, 2], "lipid_name": ["PC 34:1", "DG 34:1"], "lipid_class": ["PC", "DG"],
        "confidence": ["Gold", "Silver"], "quant_route": ["Q1", "Q3"],
        "is_compound": ["PC_d7", "TG_d7"], "RIE_pred": [np.nan, np.nan]})
    area = pd.DataFrame({"s1": [200.0, 100.0]}, index=[1, 2])
    is_area = pd.DataFrame({"s1": [100.0, 50.0]}, index=["PC_d7", "TG_d7"])
    samples = pd.DataFrame({"file_id": ["s1"], "dry_mass_g": [2.0], "spike_ng": [10.0]})
    cal = pd.DataFrame({"compound": ["PC_d7", "TG_d7"], "cal_min_area": [10.0, 100.0],
                        "cal_max_area": [1000.0, 1000.0]})
    out = estimate(routes, area, is_area, samples, standards, rie, cal, polarity="pos").set_index("feature_id")
    # Q1: (200/100) x 10 ng x f=1 / 2 g = 10 ng/g
    assert out.loc[1, "ng_per_g"] == pytest.approx(10.0)
    assert out.loc[1, "f_source"] == "class-matched standard (f = 1)"
    # Q3: f = RIE(TG, M+NH4) / RIE(DG, M+NH4) = 10 / 5 = 2 -> (100/50) x 10 x 2 / 2 = 20 ng/g
    assert out.loc[2, "f"] == pytest.approx(2.0)
    assert out.loc[2, "ng_per_g"] == pytest.approx(20.0)
    assert out.loc[2, "is_range_flag"] == "below"


def test_ammonium_responses_are_bridged_to_the_lpe_scale_through_ps():
    rie = pd.DataFrame({"class": ["PS", "PS", "TG", "LPE"], "adduct": ["M+H", "M+NH4", "M+NH4", "M+H"],
                        "RIE_LPE": [7.0, 1.0, 13.0, 1.0]})
    r = class_rie(rie)
    assert r[("TG", "M+NH4")] == pytest.approx(13.0 * 7.0)
    assert r[("PS", "M+NH4")] == pytest.approx(7.0)
    assert r[("LPE", "M+H")] == 1.0


def test_formate_responses_are_bridged_through_hexcer():
    rie = pd.DataFrame({"class": ["HexCer", "HexCer", "PC"], "adduct": ["M-H", "M+HCOO", "M+HCOO"],
                        "RIE_LPE": [0.25, 1.0, 2.0]})
    assert class_rie(rie)[("PC", "M+HCOO")] == pytest.approx(0.5)


def test_f_is_capped_and_flagged(standards):
    rie = pd.DataFrame({"class": ["TG", "DG", "PS", "PS"], "adduct": ["M+NH4", "M+NH4", "M+H", "M+NH4"],
                        "RIE_LPE": [100.0, 1.0, 1.0, 1.0]})
    routes = pd.DataFrame({"feature_id": [2], "lipid_name": ["DG 34:1"], "lipid_class": ["DG"],
                           "confidence": ["Silver"], "quant_route": ["Q3"], "is_compound": ["TG_d7"],
                           "RIE_pred": [np.nan]})
    out = estimate(routes, pd.DataFrame({"s1": [10.0]}, index=[2]), pd.DataFrame({"s1": [10.0]}, index=["TG_d7"]),
                   pd.DataFrame({"file_id": ["s1"], "dry_mass_g": [1.0], "spike_ng": [10.0]}), standards, rie,
                   pd.DataFrame(columns=["compound", "cal_min_area", "cal_max_area"]), polarity="pos")
    assert out.loc[0, "f"] == 10.0 and bool(out.loc[0, "f_capped"])
    assert out.loc[0, "f_uncapped"] == pytest.approx(100.0)
