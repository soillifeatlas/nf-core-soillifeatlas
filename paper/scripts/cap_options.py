#!/usr/bin/env python3
"""Review item 1: what-if for the Q3 response-factor cap (positive mode). Changes nothing in stage 5.

Option A (fix): score each analyte at its class's usual positive-mode adduct
  ([M+NH4]+ for TG, DG, MG, CE, PI, PG; [M+H]+ otherwise) instead of the standard's adduct,
  with the same PS bridge and the same 0.1-10 cap.
Option B (demote): Q3 features of the affected classes (TG, DG, MG, PI, PG) become Q4
  (relative only, excluded from mass sums).
Amounts scale linearly with f, so new ng/g = current ng/g x f_new / f_current.

Writes diagnosis_f_by_class.tsv, class_totals_by_option.tsv, treatment_totals_by_option.tsv,
plot_totals_by_option.tsv.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from soillifeatlas import quantify as qt  # noqa: E402
import provenance  # noqa: E402

AFFECTED = ["TG", "DG", "MG", "PI", "PG"]
NH4_POS = {"TG", "DG", "MG", "CE", "PI", "PG"}


def f_option_a(row, std, rie_m):
    if row["quant_route"] != "Q3":
        return row["f"]
    k = row["is_compound"]
    s_cls = qt._cls(std.at[k, "lipid_class"])
    s_add = qt._adduct_key(std.at[k, "adduct"])
    a_cls = qt._cls(row["lipid_class"])
    a_add = "M+NH4" if a_cls in NH4_POS else "M+H"
    r_s, r_a = rie_m.get((s_cls, s_add)), rie_m.get((a_cls, a_add))
    if r_s is None or r_a is None:
        return row["f"]
    return float(np.clip(r_s / r_a, qt.F_MIN, qt.F_MAX))


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--quant", required=True, help="paper/04_quant/pos")
    ap.add_argument("--rie", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args(argv)
    q, out = Path(a.quant), Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    est_p, per_p, std_p = q / "estimates_long.tsv", q / "per_lipid.tsv", q / "standards_used.tsv"
    est, per = pd.read_csv(est_p, sep="\t"), pd.read_csv(per_p, sep="\t")
    std = pd.read_csv(std_p, sep="\t").set_index("compound")
    rie_m = qt.class_rie(pd.read_csv(a.rie))

    # diagnosis per class
    rows = []
    for c in AFFECTED + ["PE", "PC", "SM", "LPC", "CER"]:
        rows.append({"class": c, "RIE_M+H": rie_m.get((c, "M+H")), "RIE_M+NH4_table_PS=1": None,
                     "RIE_M+NH4_bridged": rie_m.get((c, "M+NH4"))})
    raw = pd.read_csv(a.rie)
    raw["c"] = raw["class"].str.upper()
    rawm = raw.groupby(["c", "adduct"])["RIE_LPE"].mean()
    for r in rows:
        r["RIE_M+NH4_table_PS=1"] = rawm.get((r["class"], "M+NH4"))
    d = pd.DataFrame(rows)
    d["bridge_factor_PS"] = rawm[("PS", "M+H")] / rawm[("PS", "M+NH4")]
    q3 = per[per["quant_route"] == "Q3"]
    g = q3.groupby(q3["lipid_class"].str.upper()).agg(n_q3=("f", "size"), pct_capped=("f_capped", "mean"),
                                                       median_f_uncapped=("f_uncapped", "median"))
    d = d.merge(g, left_on="class", right_index=True, how="left")
    d["pct_capped"] = (d["pct_capped"] * 100).round(1)
    d.to_csv(out / "diagnosis_f_by_class.tsv", sep="\t", index=False)

    # options
    fa = per.apply(lambda r: f_option_a(r, std, rie_m), axis=1)
    per = per.assign(f_A=fa)
    e = est.merge(per[["feature_id", "f_A", "lipid_id"]], on="feature_id")
    e["ng_current"] = e["ng_per_g"]
    e["ng_A"] = e["ng_per_g"] * e["f_A"] / e["f"]
    e["ng_B"] = np.where((e["quant_route"] == "Q3") & e["lipid_class"].str.upper().isin(AFFECTED), np.nan,
                         e["ng_per_g"])
    keep = e.groupby(["lipid_id", "feature_id"])["area"].median().reset_index() \
            .sort_values("area", ascending=False).drop_duplicates("lipid_id")["feature_id"]
    e = e[e["feature_id"].isin(keep)]
    cols = ["ng_current", "ng_A", "ng_B"]
    by_plot_class = e.groupby(["plot", "treatment", "lipid_class"])[cols].sum(min_count=1).reset_index()
    cls = by_plot_class.groupby("lipid_class")[cols].mean() / 1000
    cls.columns = ["ug_g_current", "ug_g_A_fix", "ug_g_B_demote"]
    cls.sort_values("ug_g_current", ascending=False).round(3).to_csv(out / "class_totals_by_option.tsv", sep="\t")
    plot = by_plot_class.groupby(["plot", "treatment"])[cols].sum(min_count=1).reset_index()
    for c in cols:
        plot[c] = plot[c] / 1000
    plot.round(2).to_csv(out / "plot_totals_by_option.tsv", sep="\t", index=False)
    tr = plot.groupby("treatment")[cols].agg(["mean", "std"]).round(1)
    tr.to_csv(out / "treatment_totals_by_option.tsv", sep="\t")
    per[["feature_id", "lipid_name", "lipid_class", "quant_route", "is_compound", "f", "f_uncapped", "f_A"]].to_csv(
        out / "per_lipid_f_option_A.tsv", sep="\t", index=False)
    provenance.write(out, {"estimates_long": est_p, "per_lipid": per_p, "standards_used": std_p, "rie": a.rie},
                     {"affected_classes": AFFECTED, "option_A_nh4_classes": sorted(NH4_POS),
                      "cap": [qt.F_MIN, qt.F_MAX]})


if __name__ == "__main__":
    main()
