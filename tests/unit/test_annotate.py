"""Unit tests for stage 4 annotation merging (soillifeatlas.annotate)."""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from soillifeatlas.annotate import (
    check_composition,
    check_rt,
    from_library,
    from_lipidsearch,
    from_mzmine,
    from_sirius,
    ion_mz,
    match_by_mz_rt,
    merge,
    parse_lipid,
    summarise,
)


@pytest.mark.parametrize("text, cls, level, name, summ", [
    ("PC 16:0_18:1", "PC", "species", "PC 16:0_18:1", "34:1"),
    ("PC 34:1", "PC", "sum", "PC 34:1", "34:1"),
    ("TG(16:0_16:0_16:1)", "TG", "species", "TG 16:0_16:0_16:1", "48:1"),
    ("LPC 18:1", "LPC", "species", "LPC 18:1", "18:1"),
    ("Cer d18:1/15:0", "Cer", "species", "Cer d18:1_15:0", "d33:1"),
    ("PE O-16:0_18:1", "PE", "species", "PE O-16:0_18:1", "O-34:1"),
    ("DG(O-25:2)", "DG", "sum", "DG O-25:2", "O-25:2"),
    ("PC(16:0/18:1); [M+H]+ C42H83N1O8P1", "PC", "species", "PC 16:0_18:1", "34:1"),
    ("Massbank:LQB00208 PE 36:2", "PE", "sum", "PE 36:2", "36:2"),
    ("Spectral Match to Lyso-PC(16:0) from NIST14", "LPC", "species", "LPC 16:0", "16:0"),
    ("PC 18:2", "LPC", "species", "LPC 18:2", "18:2"),  # one short chain = lyso form
    ("PE(18:1/0:0)", "LPE", "species", "LPE 18:1", "18:1"),
])
def test_parse_lipid(text, cls, level, name, summ):
    lip = parse_lipid(text)
    assert lip is not None
    assert lip.lipid_class == cls
    assert lip.level == level
    assert lip.name() == name
    assert lip.sum_composition == summ


@pytest.mark.parametrize("text", ["DEHP CollisionEnergy:205060", "Spectral Match to 1-Palmitoyl-2-hydroxy from NIST14", "", None])
def test_parse_rejects_non_lipids(text):
    assert parse_lipid(text) is None


def test_ion_mz():
    assert ion_mz(752.6061, "M+H") == pytest.approx(753.6134, abs=1e-4)
    assert ion_mz(752.6061, "[M+NH4]+") == pytest.approx(770.6399, abs=1e-4)
    assert ion_mz(752.6061, "M-H") == pytest.approx(751.5988, abs=1e-4)
    assert ion_mz(100.0, "M+Xyz") is None


def _features(n=40):
    rng = np.random.default_rng(1)
    mz = np.sort(rng.uniform(400, 900, n))
    rt = rng.uniform(2, 25, n)
    return pd.DataFrame({"id": np.arange(1, n + 1), "mz": mz, "rt": rt})


def test_match_removes_rt_offset():
    f = _features()
    ext = pd.DataFrame({"mz": f["mz"] * (1 + 2e-6), "rt": f["rt"] - 0.35})  # 0.35 min earlier in the other run
    out = match_by_mz_rt(f, ext, ppm=5, rt_window=0.1)
    assert out.attrs["rt_offset"] == pytest.approx(0.35, abs=0.01)
    assert (out["feature_id"].to_numpy() == f["id"].to_numpy()).all()


def test_lipidsearch_grades():
    f = pd.DataFrame({"id": [1, 2, 3], "mz": [753.6134, 776.5538, 851.7598], "rt": [15.7, 16.1, 24.0]})
    ls = pd.DataFrame({
        "Rej": [False, False, False],
        "LipidMolec": ["PC(15:0_18:1)", "PE(18:1_20:4)", "TG(16:0_18:1_17:0)"],
        "CalcMass": [752.6061, 775.5465, 850.7525],
        "BaseRt": [15.72, 16.12, 24.01],
        "MainIon": ["M+H", "M+H", "M+H"],
        "TotalGrade": ["A", "C", "B"],
    })
    c = from_lipidsearch(ls, f, ppm=10, rt_window=0.2).set_index("feature_id")
    assert c.loc[1, "rank"] == 1 and c.loc[1, "level"] == "species"
    assert c.loc[2, "rank"] == 3 and c.loc[2, "level"] == "sum"  # grade C: sum level at most
    assert c.loc[3, "rank"] == 1


def test_library_thresholds_and_same_run():
    f = _features()
    lib = pd.DataFrame({
        "Compound_Name": ["PC(16:0/18:1); [M+H]+ C42H83N1O8P1", "PE 34:1", "DEHP"],
        "MQScore": [0.9, 0.5, 0.95],
        "SharedPeaks": [8, 10, 9],
        "#Scan#": [5, 6, 7],
    })
    c = from_library(lib, f, 5, 0.2, min_score=0.7, min_peaks=4, same_run=True)
    assert list(c["feature_id"]) == [5]
    assert c.iloc[0]["rank"] == 2


def test_sirius_class_only():
    f = _features()
    s = pd.DataFrame({"featureId": [3, 4], "ClassyFire#most specific class": ["Glycerophosphocholines", "Benzenoids"],
                      "ClassyFire#most specific class Probability": [0.95, 0.99]})
    c = from_sirius(s, f, 5, 0.2, same_run=True)
    assert list(c["feature_id"]) == [3]
    assert c.iloc[0]["level"] == "class" and c.iloc[0]["rank"] == 4 and c.iloc[0]["lipid_class"] == "PC"


def test_merge_confidence_rules():
    f = pd.DataFrame({"id": [1, 2, 3, 4], "mz": [1, 2, 3, 4.0], "rt": [1, 2, 3, 4.0]})
    full = pd.DataFrame({"id": [1, 2], "lipid_annotations:lipid_annotations": ["PC 16:0_18:1", "PE 34:1"],
                         "lipid_annotations:ion_adduct": ["[M+H]+", "[M+H]+"],
                         "lipid_annotations:explained_intensity_percent": [0.8, 0.7]})
    mz = from_mzmine(full)
    lib = from_library(pd.DataFrame({"Compound_Name": ["PE(16:0/18:1)", "TG(16:0/18:1/18:1)"], "MQScore": [0.9, 0.9],
                                     "SharedPeaks": [8, 8], "#Scan#": [2, 3]}), f, 5, 0.2, same_run=True)
    sir = from_sirius(pd.DataFrame({"featureId": [4], "ClassyFire#subclass": ["Triradylcglycerols"],
                                    "ClassyFire#subclass Probability": [0.9]}), f, 5, 0.2, same_run=True)
    cands = pd.concat([mz, lib, sir], ignore_index=True)
    m = merge(cands, f, is_classes={"PC"}).set_index("feature_id")
    assert m.loc[1, "confidence"] == "Gold" and m.loc[1, "quant_route"] == "Q1"  # species from lipid rules, PC standard present
    assert m.loc[2, "confidence"] == "Silver" and m.loc[2, "lipid_name"] == "PE 34:1"  # library species, rules sum
    assert m.loc[3, "confidence"] == "Silver" and m.loc[3, "quant_route"] == "Q2"  # library species alone -> sum
    assert m.loc[4, "confidence"] == "Bronze" and m.loc[4, "quant_route"] == "Q3"
    s = summarise(m.reset_index(), cands, f)
    assert s["confidence"] == {"Gold": 1, "Silver": 2, "Bronze": 1}


def test_two_sources_agreeing_on_species_make_gold():
    f = pd.DataFrame({"id": [1], "mz": [760.585], "rt": [16.0]})
    lib = from_library(pd.DataFrame({"Compound_Name": ["PC(16:0/18:1)"], "MQScore": [0.9], "SharedPeaks": [8],
                                     "#Scan#": [1]}), f, 5, 0.2, same_run=True)
    ls = from_lipidsearch(pd.DataFrame({"LipidMolec": ["PC(16:0_18:1)"], "CalcMass": [759.5778], "BaseRt": [16.0],
                                        "MainIon": ["M+H"], "TotalGrade": ["C"]}), f, 5, 0.2)
    ls["level"] = "species"  # a species-level claim from a second source
    ls["lipid_name"] = "PC 16:0_18:1"
    m = merge(pd.concat([lib, ls], ignore_index=True), f)
    assert m.iloc[0]["confidence"] == "Gold"


def test_lyso_and_family_agree():
    f = pd.DataFrame({"id": [1], "mz": [496.34], "rt": [3.0]})
    lib = from_library(pd.DataFrame({"Compound_Name": ["PC 16:0"], "MQScore": [0.9], "SharedPeaks": [8],
                                     "#Scan#": [1]}), f, 5, 0.2, same_run=True)
    sir = from_sirius(pd.DataFrame({"featureId": [1], "ClassyFire#most specific class": ["Glycerophosphocholines"],
                                    "ClassyFire#most specific class Probability": [0.9]}), f, 5, 0.2, same_run=True)
    m = merge(pd.concat([lib, sir], ignore_index=True), f)
    assert m.iloc[0]["lipid_name"] == "LPC 16:0" and bool(m.iloc[0]["class_agreement"])


def test_check_rt_flags_off_trend():
    # TG retention rises 0.5 min per carbon and falls 0.4 min per double bond
    comps = [(48, 1), (50, 1), (50, 2), (52, 2), (52, 3), (54, 3), (54, 4), (56, 2)]
    rows = [{"feature_id": i, "rt": 0.5 * c - 0.4 * d - 8 + 0.05 * (i % 2), "lipid_name": f"TG {c}:{d}",
             "lipid_class": "TG", "confidence": "Gold"} for i, (c, d) in enumerate(comps)]
    rows += [{"feature_id": 90, "rt": 1.4, "lipid_name": "TG 50:2", "lipid_class": "TG", "confidence": "Silver"},
             {"feature_id": 91, "rt": 16.2, "lipid_name": "TG 50:2", "lipid_class": "TG", "confidence": "Silver"},
             {"feature_id": 92, "rt": 1.5, "lipid_name": "TG 26:4", "lipid_class": "TG", "confidence": "Silver"},
             {"feature_id": 93, "rt": 5.0, "lipid_name": "PI", "lipid_class": "PI", "confidence": "Bronze"}]
    out = check_rt(pd.DataFrame(rows)).set_index("feature_id")
    assert out.loc[90, "rt_check"] == "outlier"
    assert out.loc[91, "rt_check"] == "ok"
    assert out.loc[92, "rt_check"] == "untested"  # far outside the anchors' carbon range
    assert out.loc[93, "rt_check"] == "untested"
    assert (out.loc[range(8), "rt_check"] == "ok").all()


@pytest.mark.parametrize("name, cls, verdict", [
    ("TG 16:0_18:1_18:1", "TG", "ok"),
    ("TG 26:4", "TG", "implausible"),        # under 10 carbons per chain
    ("DG 7:0_6:0", "DG", "implausible"),     # short named chains
    ("DG P-4:0_17:0", "DG", "implausible"),
    ("TG 54:19", "TG", "implausible"),       # over 6 double bonds per chain
    ("PC 36:2", "PC", "ok"),
    ("LPC 16:0", "LPC", "ok"),
    ("Cer d18:1_16:0", "Cer", "ok"),
    ("PC", "PC", "untested"),
    ("AEG 34:1", "AEG", "untested"),         # chain count unknown
    ("TG 9:0_16:0_18:1+O", "TG", "ok"),      # truncated oxidised chain
    ("TG 2:0_15:0_17:1+OO", "TG", "implausible"),
])
def test_check_composition(name, cls, verdict):
    assert check_composition(name, cls) == verdict


def test_check_rt_fits_ether_apart_and_ignores_implausible_anchors():
    rows = []
    for i, (c, d) in enumerate([(32, 0), (34, 1), (34, 2), (36, 2), (36, 3), (38, 4), (40, 6)]):
        rows.append({"feature_id": i, "rt": 0.6 * c - 0.5 * d - 5, "lipid_name": f"PC {c}:{d}",
                     "lipid_class": "PC", "confidence": "Gold"})
    # an implausible Gold anchor at an absurd time must not widen the fit
    rows.append({"feature_id": 50, "rt": 1.0, "lipid_name": "PC 9:0_9:0", "lipid_class": "PC", "confidence": "Gold"})
    # ether species elute later; with too few ether anchors they stay untested, not outliers
    rows.append({"feature_id": 60, "rt": 17.0, "lipid_name": "PC O-34:1", "lipid_class": "PC", "confidence": "Silver"})
    rows.append({"feature_id": 61, "rt": 9.0, "lipid_name": "PC 34:1", "lipid_class": "PC", "confidence": "Silver"})
    out = check_rt(pd.DataFrame(rows)).set_index("feature_id")
    assert out.loc[50, "composition_check"] == "implausible"
    assert (out.loc[range(7), "rt_check"] == "ok").all()
    assert out.loc[60, "rt_check"] == "untested"
    assert out.loc[61, "rt_check"] == "outlier"  # expected 14.9 min
