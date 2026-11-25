"""Stage 5: semi-quantitative lipid amounts per gram dry soil (Q1-Q4).

For analyte j in sample i quantified against internal standard k
(manuscript section 2.7):

    c(i,j) = A(i,j) / A(i,k) x m_spiked(i,k) x f(j,k) / M(i)

with f = standard response / analyte response on the same amount basis.
No separate recovery factor is applied: losses shared by analyte and
standard cancel in the ratio.

Routes follow the response evidence available for each feature (decision
2026-11-25, option A of paper/04_quant/DESIGN_NOTE.md):

- Q1  a standard of the same lipid class is in the run       -> that standard, f = 1
- Q2  no class-matched standard, but a lipid-specific predicted
      response exists (Gold/Silver name with a prediction)    -> nearest standard,
                                                                 f = RIE(standard class) / RIE_pred
- Q3  class known, no Q1/Q2                                   -> nearest standard in RT,
                                                                 f = RIE(standard class) / RIE(analyte class)
- Q4  no lipid class                                          -> no mass; excluded from sums

Annotation confidence stays a separate column and does not decide the route.
RIE values are class means per adduct from the response table (Samrat et al.
2025, Table S10), with [M+NH4]+ and [M+HCOO]- rows bridged onto the LPE scale
(see BRIDGES). f is limited to [F_MIN, F_MAX]; ``f_capped`` marks limited rows
and ``f_uncapped`` keeps the raw ratio. The standard's measured adduct is used for both classes; if
the analyte class has no value for that adduct, its [M+H]+ / [M-H]- value is
used (``f_source`` says so); with no value at all f = 1 and ``f_default`` is set.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from .annotate import class_family

ROUTE_ORDER = ["Q1", "Q2", "Q3", "Q4"]


def _adduct_key(adduct: str) -> str:
    """'[M+NH4]+' -> 'M+NH4' (the response table's spelling)."""
    a = str(adduct).strip()
    if a.startswith("["):
        a = a[1:a.rfind("]")]
    return a


def _cls(c) -> str | None:
    if c is None or (isinstance(c, float) and np.isnan(c)) or str(c).strip() == "":
        return None
    return class_family(str(c)).upper()


# Table S10 normalises [M+NH4]+ to PS = 1 and [M+HCOO]- to HexCer = 1, not to LPE.
# Bring them onto the LPE scale of the protonated / deprotonated rows through the
# class measured with both adducts (decision 2026-11-25).
BRIDGES = {"M+NH4": ("PS", "M+H"), "M+HCOO": ("HEXCER", "M-H")}
F_MIN, F_MAX = 0.1, 10.0  # decision 2026-11-25: cap the response correction, flag capped rows


def class_rie(rie: pd.DataFrame) -> dict[tuple[str, str], float]:
    """Mean RIE per (class, adduct) on one LPE scale; class names upper-cased for matching."""
    t = rie.assign(_c=rie["class"].map(lambda c: str(c).upper()), _a=rie["adduct"].map(_adduct_key))
    m = t.groupby(["_c", "_a"])["RIE_LPE"].mean().to_dict()
    for adduct, (anchor, base) in BRIDGES.items():
        if (anchor, adduct) in m and (anchor, base) in m:
            k = m[(anchor, base)] / m[(anchor, adduct)]
            for key in [key for key in m if key[1] == adduct]:
                m[key] = m[key] * k
    return m


def assign_routes(ann: pd.DataFrame, standards: pd.DataFrame, predicted: pd.DataFrame | None) -> pd.DataFrame:
    """One row per annotated feature: route, standard used and lipid-specific RIE if any.

    ``standards``: compound, lipid_class, rt (median RT in this run), adduct.
    ``predicted``: optional lipid_name, RIE_pred.
    """
    std = standards.assign(_c=standards["lipid_class"].map(_cls)).dropna(subset=["rt"])
    by_class = dict(zip(std["_c"], std["compound"]))
    pred = {} if predicted is None else dict(zip(predicted["lipid_name"], predicted["RIE_pred"]))
    rows = []
    for _, r in ann.iterrows():
        c = _cls(r.get("lipid_class"))
        named = r.get("confidence") in ("Gold", "Silver")
        rie_pred = pred.get(r["lipid_name"], np.nan) if named else np.nan
        if c is None:
            route, k = "Q4", None
        elif c in by_class:
            route, k = "Q1", by_class[c]
        else:
            route = "Q2" if np.isfinite(rie_pred) else "Q3"
            k = std.iloc[(std["rt"] - r["rt"]).abs().argmin()]["compound"] if len(std) else None
        rows.append({"feature_id": r["feature_id"], "lipid_name": r["lipid_name"], "lipid_class": r.get("lipid_class"),
                     "confidence": r.get("confidence"), "quant_route": route, "is_compound": k,
                     "RIE_pred": rie_pred if route == "Q2" else np.nan})
    return pd.DataFrame(rows)


def range_flag(is_area: float, cal_min: float, cal_max: float) -> str:
    """Where the standard's area in a sample lies relative to its neat calibration areas."""
    if not np.isfinite(is_area) or is_area <= 0:
        return "is_not_detected"
    if np.isfinite(cal_min) and is_area < cal_min:
        return "below"
    if np.isfinite(cal_max) and is_area > cal_max:
        return "above"
    return "within" if np.isfinite(cal_min) else "no_calibration"


def _f(route, analyte_cls, std_cls, std_adduct, rie_pred, rie, polarity):
    if route == "Q1":
        return 1.0, "class-matched standard (f = 1)", False
    a = _adduct_key(std_adduct)
    r_std = rie.get((std_cls, a))
    if route == "Q2":
        if r_std is None:
            return 1.0, "no RIE for standard class; f = 1", True
        return r_std / rie_pred, f"RIE({std_cls},{a}) / predicted RIE", False
    r_an, src = rie.get((analyte_cls, a)), f"RIE({std_cls},{a}) / RIE({analyte_cls},{a})"
    if r_an is None:
        base = "M+H" if polarity == "pos" else "M-H"
        r_an, src = rie.get((analyte_cls, base)), f"RIE({std_cls},{a}) / RIE({analyte_cls},{base}) [adduct fallback]"
    if r_std is None or r_an is None:
        return 1.0, "class not in response table; f = 1", True
    return r_std / r_an, src, False


def estimate(routes: pd.DataFrame, area: pd.DataFrame, is_area: pd.DataFrame, samples: pd.DataFrame,
             standards: pd.DataFrame, rie: pd.DataFrame, calibration: pd.DataFrame,
             polarity: str) -> pd.DataFrame:
    """Long table: one row per feature x sample with ng per g dry soil and its flags.

    ``area``: features x sample files (peak areas). ``is_area``: standard compound x
    sample files (NaN where not detected). ``samples``: file_id, dry_mass_g, spike_ng.
    ``calibration``: compound, cal_min_area, cal_max_area from the neat calibration.
    """
    rie_m = class_rie(rie)
    std = standards.set_index("compound")
    cal = calibration.set_index("compound")
    out = []
    for _, r in routes.iterrows():
        if r["quant_route"] == "Q4" or r["is_compound"] is None or r["feature_id"] not in area.index:
            continue
        k = r["is_compound"]
        f_raw, f_src, f_def = _f(r["quant_route"], _cls(r["lipid_class"]), _cls(std.at[k, "lipid_class"]),
                                 std.at[k, "adduct"], r.get("RIE_pred", np.nan), rie_m, polarity)
        f = float(np.clip(f_raw, F_MIN, F_MAX))
        f_capped = f != f_raw
        for _, s in samples.iterrows():
            fid = s["file_id"]
            a_j = area.at[r["feature_id"], fid] if fid in area.columns else np.nan
            a_k = is_area.at[k, fid] if (k in is_area.index and fid in is_area.columns) else np.nan
            ok = np.isfinite(a_j) and a_j > 0 and np.isfinite(a_k) and a_k > 0
            ng_g = a_j / a_k * s["spike_ng"] * f / s["dry_mass_g"] if ok else np.nan
            out.append({
                "feature_id": r["feature_id"], "lipid_name": r["lipid_name"], "lipid_class": r["lipid_class"],
                "confidence": r["confidence"], "quant_route": r["quant_route"], "file_id": fid,
                "is_compound": k, "f": f, "f_uncapped": f_raw, "f_capped": f_capped,
                "f_source": f_src, "f_default": f_def,
                "area": a_j, "is_area": a_k, "ng_per_g": ng_g,
                "is_range_flag": range_flag(a_k, cal["cal_min_area"].get(k, np.nan),
                                            cal["cal_max_area"].get(k, np.nan)),
            })
    return pd.DataFrame(out)
