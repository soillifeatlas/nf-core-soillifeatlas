"""Figure 5 - response, recovery, routes and amounts."""
from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from common import (DOUBLE, MM, POL, POL_LABEL, Panels, fig_letter, gates, input_path, letter,
                    pending, read_table, save, set_panel, show_image)

POLS = ("pos", "neg")
ADDUCT = {  # adduct -> (polarity, marker, label)
    "M+H": ("pos", "o", r"[M+H]$^{+}$"), "M+NH4": ("pos", "s", r"[M+NH$_4$]$^{+}$"),
    "M-H": ("neg", "o", r"[M−H]$^{-}$"), "M+HCOO": ("neg", "s", r"[M+HCOO]$^{-}$"),
}
SCHEM = "figures/schematics/source/Fig5_quantification_2026-11-24.png"
CROP_5D = (185, 1363, 3664, 2087)


def panel_A(ax, tag="A"):
    set_panel("5A")
    d = read_table("03_response/response_by_class_adduct.tsv")
    unknown = set(d.adduct) - set(ADDUCT)
    if unknown:
        raise ValueError(f"unmapped adducts: {unknown}")
    order = d.groupby("class").logRIE_LPE.median().sort_values(ascending=False).index
    xpos = {c: i for i, c in enumerate(order)}
    rng = np.random.default_rng(0)  # horizontal jitter only, for legibility
    for add, (pol, mk, lab) in ADDUCT.items():
        s = d[d.adduct == add]
        x = s["class"].map(xpos) + rng.uniform(-0.18, 0.18, len(s))
        ax.scatter(x, s.RIE_LPE, marker=mk, s=10, lw=0.5,
                   facecolor=POL[pol] if mk == "o" else "white", edgecolor=POL[pol],
                   label=f"{lab} (n = {len(s)})", zorder=3)
    ax.set_yscale("log")
    ax.axhline(1, color="#888888", lw=0.5, ls=":")
    ax.set_xticks(range(len(order)))
    ax.set_xticklabels(order, rotation=90)
    ax.set_xlim(-0.7, len(order) - 0.3)
    ax.set_ylabel("RIE relative to LPE (log scale)")
    ax.legend(loc="lower left", ncol=2, handletextpad=0.2, columnspacing=0.8)
    span = np.log10(d.RIE_LPE.max() / d.RIE_LPE.min())
    letter(ax, tag, x=-0.1, y=1.03,
           title=f"Response spans {span:.0f} orders of magnitude across classes and adducts "
                 f"({len(d)} standards)")


def panel_B(ax, tag="B"):
    set_panel("5B")
    p = read_table("03_response/ie_loco_predictions.tsv")
    st = read_table("03_response/ie_model_stats.tsv")
    best = st[st.best.astype(bool)].iloc[0]
    ax.scatter(p.pred_best, p.logRIE_LPE, s=9, color=POL["pos"], lw=0, alpha=0.85, zorder=3)
    lo = min(p.pred_best.min(), p.logRIE_LPE.min()) - 0.3
    hi = max(p.pred_best.max(), p.logRIE_LPE.max()) + 0.3
    ax.plot([lo, hi], [lo, hi], color="black", lw=0.6, ls="--")
    ax.set_xlim(lo, hi)
    ax.set_ylim(lo, hi)
    ax.set_aspect("equal")
    ax.set_xlabel(r"predicted log$_{10}$ RIE (leave-one-class-out)")
    ax.set_ylabel(r"measured log$_{10}$ RIE (rel. LPE)")
    ax.text(0.03, 0.97,
            f"{best.model.replace('_', ' ')}\nR² = {best.r2:.2f}\n"
            f"RMSE = {best.rmse_log10:.2f} log$_{{10}}$ (~{best.fold_error:.0f}-fold)\n"
            f"{int(best.n_standards)} POS standards, {int(best.n_classes)} classes",
            transform=ax.transAxes, va="top", fontsize=6)
    letter(ax, tag, x=-0.2, y=1.03,
           title=f"Response is not predictable across classes (R² = {best.r2:.2f})")


def panel_C(ax, tag="C"):
    set_panel("5C")
    data = {pol: read_table(f"01_qc/{pol}/is_vs_calibration.tsv") for pol in POLS}
    if not gates()["reconstitution_volume_confirmed"]:
        pending(ax, "reconstitution\nvolume not confirmed")
        letter(ax, tag, x=0.0, y=1.03, title="Apparent recovery in soil")
        return "PENDING", "gated on reconstitution_volume_confirmed"
    rows = []
    for pol, d in data.items():
        d = d.dropna(subset=["R_app_median_pct"]).copy()
        d["pol"] = pol
        rows.append(d)
    d = pd.concat(rows).reset_index(drop=True)
    for i, r in d.iterrows():
        ax.plot([i, i], [r.R_app_min_pct, r.R_app_max_pct], color=POL[r.pol], lw=0.8)
        ax.scatter(i, r.R_app_median_pct, color=POL[r.pol], s=12, zorder=3)
    ax.axhline(100, color="black", lw=0.6, ls="--")
    ax.set_xticks(range(len(d)))
    ax.set_xticklabels([f"{c.split('_')[0]} {POL_LABEL[p]}" for c, p in zip(d.compound, d.pol)],
                       rotation=90)
    ax.set_ylabel("R_app (%), median and range over soils")
    letter(ax, tag, x=-0.15, y=1.03, title="Apparent recovery in soil")
    return "REAL", ""


def panel_E(ax, tag="E"):
    set_panel("5E")
    include_q3 = bool(gates()["stage5_fix_confirmed"])
    lo, hi = np.inf, -np.inf
    lines = []
    for pol in POLS:
        d = read_table(f"05_nist/{pol}/nist_comparison.tsv")
        q1 = d[d.routes == "Q1"]
        lab = f"{POL_LABEL[pol]} Q1, n = {len(q1)}"
        if len(q1):
            lab += (f"\nmedian ratio {np.median(q1.ratio):.2f}, "
                    f"{q1.within_2fold.astype(bool).mean() * 100:.0f}% within 2-fold")
        ax.scatter(q1.medm, q1.estimate_nmol_per_ml, s=10, color=POL[pol], lw=0, zorder=3,
                   label=lab)
        sets = [q1]
        if include_q3:
            q3 = d[d.routes != "Q1"]
            ax.scatter(q3.medm, q3.estimate_nmol_per_ml, s=10, facecolor="white",
                       edgecolor=POL[pol], lw=0.6, zorder=3,
                       label=f"{POL_LABEL[pol]} Q3 (n = {len(q3)})")
            sets.append(q3)
        for s in sets:
            if len(s):
                lo = min(lo, s.medm.min(), s.estimate_nmol_per_ml.min())
                hi = max(hi, s.medm.max(), s.estimate_nmol_per_ml.max())
    lo, hi = lo / 2.5, hi * 2.5
    xx = np.array([lo, hi])
    ax.plot(xx, xx, color="black", lw=0.6)
    for f in (2, 0.5):
        ax.plot(xx, xx * f, color="black", lw=0.5, ls="--")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(lo, hi)
    ax.set_ylim(lo, hi)
    ax.set_aspect("equal")
    ax.set_xlabel(r"NIST SRM 1950 consensus (nmol mL$^{-1}$)")
    ax.set_ylabel(r"estimate (nmol mL$^{-1}$)")
    ax.legend(loc="lower right", handletextpad=0.2, fontsize=5, borderaxespad=0.2)
    ax.text(0.03, 0.97, "solid 1:1\ndashed 2-fold", transform=ax.transAxes, ha="left",
            va="top", fontsize=5.5)
    q1p = read_table("05_nist/pos/nist_comparison.tsv")
    q1p = q1p[q1p.routes == "Q1"]
    letter(ax, tag, x=-0.2, y=1.03,
           title=f"Q1 vs NIST: median ratio {np.median(q1p.ratio):.2f}, "
                 f"{q1p.within_2fold.astype(bool).mean() * 100:.0f} % within 2-fold (POS)")
    if include_q3:
        return "REAL", "Q1 and Q3"
    return "REAL", "Q1 lipids only; Q3 points gated on stage5_fix_confirmed (caption: Q3 pending)"


def panel_F(ax, top=10, tag="F"):
    set_panel("5F")
    data = {pol: read_table(f"04_quant/{pol}/class_by_sample.tsv") for pol in POLS}
    if not gates()["stage5_fix_confirmed"]:
        pending(ax, "stage-5 quantification\nfix not yet confirmed")
        letter(ax, tag, x=0.0, y=1.03, title="Class amounts per g dry soil")
        return "PENDING", "gated on stage5_fix_confirmed"
    summ = []
    for pol, d in data.items():
        g = d.groupby("lipid_class").ug_per_g.agg(["mean", "std", "count"]).reset_index()
        g["pol"] = pol
        summ.append(g)
    s = pd.concat(summ)
    order = s.groupby("lipid_class")["mean"].max().sort_values(ascending=False).head(top).index
    x = np.arange(len(order))
    for k, pol in enumerate(POLS):
        sp = s[s.pol == pol].set_index("lipid_class").reindex(order)
        ax.bar(x + (k - 0.5) * 0.4, sp["mean"], 0.4, yerr=sp["std"], color=POL[pol],
               error_kw=dict(lw=0.5, capsize=1.5), label=f"{POL_LABEL[pol]}")
    ax.set_xticks(x)
    ax.set_xticklabels(order, rotation=90)
    ax.set_yscale("log")
    ax.set_ylabel(r"µg g$^{-1}$ dry soil (mean ± SD of soils)")
    ax.legend(loc="upper right")
    letter(ax, tag, x=-0.12, y=1.03, title="Class amounts per g dry soil (common lipid set)")
    return "REAL", "common lipid set; lower bounds of class totals"




def route_tree(ax, tag):
    """Native Q1-Q4 decision tree in the shared style. Labels follow the implementation
    (lib/soillifeatlas/quantify.py, option A); Q2 is drawn but not used in this run."""
    from matplotlib.patches import FancyBboxPatch, Polygon

    from common import SOIL, letter
    set_panel("route_tree")
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 62)
    ax.set_axis_off()
    xs = [12, 37, 62, 87]
    qs = ["Annotated to\na lipid class?", "Class standard\ndetected in soil?",
          "Lipid-specific\nprediction?", "Class response\navailable?"]
    yq, yr, dw, dh = 47, 13, 11, 8
    for text, x in zip(qs, xs):
        ax.add_patch(Polygon([(x - dw, yq), (x, yq + dh), (x + dw, yq), (x, yq - dh)],
                             closed=True, facecolor="#EEF1F4", edgecolor=SOIL, lw=0.7))
        ax.text(x, yq, text, ha="center", va="center", fontsize=5.3, linespacing=1.1)
    for i, (x0, x1) in enumerate(zip(xs[:-1], xs[1:])):
        ax.annotate("", xy=(x1 - dw, yq), xytext=(x0 + dw, yq),
                    arrowprops=dict(arrowstyle="-|>", lw=0.6, color="#777777", mutation_scale=6))
        ax.text((x0 + x1) / 2, yq + 1.2, "yes" if i == 0 else "no", ha="center", va="bottom",
                fontsize=5, color="#555555")
    routes = [("Q4", "relative signal only", "#DADADA", "#333333", "no"),
              ("Q1", "class-matched standard", SOIL, "white", "yes"),
              ("Q2", "not used in this run", "#F4F4F4", "#999999", "yes"),
              ("Q3", "class-level response", "#8A96A3", "white", "yes")]
    for x, (q, desc, fc, tc, lab) in zip(xs, routes):
        ax.add_patch(FancyBboxPatch((x - 11, yr - 6), 22, 12,
                                    boxstyle="round,pad=0,rounding_size=1.5", facecolor=fc,
                                    edgecolor="#BBBBBB" if q == "Q2" else "none", lw=0.6,
                                    linestyle=(0, (2, 1.5))))
        ax.text(x, yr + 2, q, ha="center", va="center", fontsize=6.5, fontweight="bold",
                color=tc)
        ax.text(x, yr - 2.5, desc, ha="center", va="center", fontsize=5, color=tc)
        ax.annotate("", xy=(x, yr + 6.2), xytext=(x, yq - dh),
                    arrowprops=dict(arrowstyle="-|>", lw=0.6, color="#777777", mutation_scale=6))
        ax.text(x + 1, (yr + 6 + yq - dh) / 2, lab, fontsize=5, color="#555555", va="center")
    # last "no" returns to Q4 underneath the boxes
    ax.plot([xs[-1] + dw, 99.5, 99.5, xs[0]], [yq, yq, 2.5, 2.5], color="#AAAAAA", lw=0.6,
            ls=(0, (2, 1.5)))
    ax.annotate("", xy=(xs[0], yr - 6.2), xytext=(xs[0], 2.5),
                arrowprops=dict(arrowstyle="-|>", lw=0.6, color="#AAAAAA", mutation_scale=6))
    ax.text(99, yq + 1.2, "no", ha="right", va="bottom", fontsize=5, color="#555555")
    letter(ax, tag, x=0.0, y=1.0, title="Each feature gets one route; Q1–Q3 give an amount, Q4 "
           "keeps relative signal")
