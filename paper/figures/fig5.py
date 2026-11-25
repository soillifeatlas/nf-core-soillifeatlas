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


def panel_A(ax):
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
    letter(ax, "A", x=-0.1, y=1.03,
           title=f"Response by class and adduct ({len(d)} standard measurements; "
                 "Samrat et al. 2025, Table S10)")


def panel_B(ax):
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
    letter(ax, "B", x=-0.2, y=1.03, title="Measured vs predicted response")


def panel_C(ax):
    set_panel("5C")
    data = {pol: read_table(f"01_qc/{pol}/is_vs_calibration.tsv") for pol in POLS}
    if not gates()["reconstitution_volume_confirmed"]:
        pending(ax, "reconstitution\nvolume not confirmed")
        letter(ax, "C", x=0.0, y=1.03, title="Apparent recovery in soil")
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
    letter(ax, "C", x=-0.15, y=1.03, title="Apparent recovery in soil")
    return "REAL", ""


def panel_E(ax):
    set_panel("5E")
    include_q3 = bool(gates()["stage5_fix_confirmed"])
    lo, hi = np.inf, -np.inf
    lines = []
    for pol in POLS:
        d = read_table(f"05_nist/{pol}/nist_comparison.tsv")
        q1 = d[d.routes == "Q1"]
        lab = f"{POL_LABEL[pol]} Q1, n = {len(q1)}"
        if len(q1):
            lab += (f"; median ratio {np.median(q1.ratio):.2f}; "
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
    ax.legend(loc="upper left", handletextpad=0.2, fontsize=5.3, borderaxespad=0.2)
    ax.text(0.97, 0.03, "solid 1:1\ndashed 2-fold", transform=ax.transAxes, ha="right",
            va="bottom", fontsize=5.5)
    letter(ax, "E", x=-0.2, y=1.03,
           title="NIST SRM 1950: estimate vs consensus" + ("" if include_q3 else " (Q1 only)"))
    if include_q3:
        return "REAL", "Q1 and Q3"
    return "REAL", "Q1 lipids only; Q3 points gated on stage5_fix_confirmed (caption: Q3 pending)"


def panel_F(ax, top=10):
    set_panel("5F")
    data = {pol: read_table(f"04_quant/{pol}/class_by_sample.tsv") for pol in POLS}
    if not gates()["stage5_fix_confirmed"]:
        pending(ax, "stage-5 quantification\nfix not yet confirmed")
        letter(ax, "F", x=0.0, y=1.03, title="Class amounts per g dry soil")
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
    letter(ax, "F", x=-0.12, y=1.03, title="Class amounts per g dry soil (common lipid set)")
    return "REAL", "common lipid set; lower bounds of class totals"


def build():
    P = Panels(5)
    fig = plt.figure(figsize=(DOUBLE, 190 * MM))
    outer = fig.add_gridspec(3, 1, height_ratios=[1.05, 0.62, 1.0], hspace=0.42,
                             left=0.08, right=0.98, top=0.95, bottom=0.09)
    r1 = outer[0].subgridspec(1, 2, width_ratios=[1.55, 1], wspace=0.28)
    panel_A(fig.add_subplot(r1[0]))
    P.add("A", "REAL", "axis: RIE relative to LPE (not response per pmol)")
    panel_B(fig.add_subplot(r1[1]))
    P.add("B", "REAL", "R², RMSE from ie_model_stats.tsv (best model)")
    set_panel("5D")
    axD = fig.add_subplot(outer[1])
    show_image(axD, input_path(SCHEM), CROP_5D)
    pos = axD.get_position()
    fig_letter(fig, "D", 0.01, pos.y1 + 0.035, "Quantification route per feature")
    P.add("D", "SCHEMATIC", "crop of Fig5_quantification_2026-11-24.png, unchanged; Q2 drawn "
          "but not used in this run (0 features)")
    r3 = outer[2].subgridspec(1, 3, width_ratios=[0.8, 1, 1.1], wspace=0.4)
    st, note = panel_C(fig.add_subplot(r3[0]))
    P.add("C", st, note)
    st, note = panel_E(fig.add_subplot(r3[1]))
    P.add("E", st, note)
    st, note = panel_F(fig.add_subplot(r3[2]))
    P.add("F", st, note)
    return save(fig, 5, "quantification"), P.rows
