#!/usr/bin/env python
"""Draw every manuscript figure from the outputs of analysis.py.

    python paper/figures/plots.py        (normally run through build_all.py)

Reads only figures/derived/ (tables and numbers.tsv written by analysis.py), plus the
report HTML it renders and the schematic source image. No number is computed here: every value
in a panel or title comes from derived/. Main text in Results order: Fig 1 Atlas (unchanged
here), Fig 2 workflow, Fig 3 what the Atlas adds, Fig 4 single features, Fig 5 ClimGrass
demonstration; supplement S1-S5. A figure that needs an unconfirmed input (gates.yaml) is not
drawn; its STATUS row says NOT BUILT.
"""
from __future__ import annotations

import os
import shutil
import subprocess

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Polygon, Rectangle

import common
from common import (ATLAS, CONF, DOUBLE, FIG_DIR, GROUP_COLOR, GROUP_ORDER, MM, NAMED, NOMATCH,
                    OI, POL, POL_LABEL, ROLE, ROLE_LABEL, SOIL, Panels, fig_letter, finding,
                    gates, input_path, letter, read_table, save, set_panel, show_image)

POLS = ("pos", "neg")


# =========================================================================== reading derived/
def D(name: str) -> pd.DataFrame:
    return read_table(f"figures/derived/{name}.tsv")


_NUM: dict[str, str] = {}


def N(key: str):
    """A number from derived/numbers.tsv, recorded against the current panel."""
    if not _NUM:
        t = pd.read_csv(common.PAPER / "figures/derived/numbers.tsv", sep="\t", dtype=str)
        _NUM.update(dict(zip(t.key, t.value)))
    if key not in _NUM:
        raise KeyError(f"number not made by analysis.py: {key}")
    common._READS.setdefault(common._current_panel[0], set()).add(
        f"figures/derived/numbers.tsv#{key}")
    v = _NUM[key]
    try:
        f = float(v)
        return int(f) if f.is_integer() and "." not in v else f
    except ValueError:
        return v


def _phyla(n: int) -> str:
    return f"{n} phylum" if n == 1 else f"{n} phyla"


def _not_built(figure, name, role, why, panels):
    P = Panels(figure, name, role)
    for p in panels:
        P.add(p, "NOT BUILT", why)
    return [], P.rows


# =========================================================================== quality panels
def _spans(pol):
    b = D("qc_blocks")
    b = b[b.pol == pol]
    return {int(r.block): (int(r.lo), int(r.hi)) for _, r in b.iterrows()}, b


def _shade(ax, spans):
    for b, (lo, hi) in spans.items():
        if b % 2 == 1:
            ax.axvspan(lo - 0.5, hi + 0.5, color="#EFEFEF", zorder=0, lw=0)


def qc_sequence(fig, gs, tag):
    set_panel("qc_sequence")
    seq = D("qc_sequence")
    axes = [fig.add_subplot(gs[i]) for i in range(2)]
    for ax, pol in zip(axes, POLS):
        inj = seq[seq.pol == pol]
        for _, r in inj.iterrows():
            ax.add_patch(Rectangle((r.injection_order - 0.45, 0), 0.9, 1,
                                   facecolor=ROLE[r.role], edgecolor="#333333", linewidth=0.3))
        spans, b = _spans(pol)
        for _, r in b.iterrows():
            ax.text((r.lo + r.hi) / 2, 1.12, r.label, ha="center", va="bottom", fontsize=6)
            if r.lo > 1:
                ax.axvline(r.lo - 0.5, color="black", lw=0.8)
        ax.set_xlim(0.4, inj.injection_order.max() + 0.6)
        ax.set_ylim(0, 1)
        ax.set_yticks([0.5])
        ax.set_yticklabels([f"{POL_LABEL[pol]} (n = {len(inj)})"])
        ax.tick_params(axis="y", length=0)
        ax.spines["left"].set_visible(False)
    axes[1].set_xlabel("injection order", labelpad=1)
    letter(axes[0], tag, x=-0.07, y=1.55,
           title=f"Calibration run first; soils follow {N('qc_block_gap_h_pos'):.0f} h later, "
                 f"bracketed by {N('qc_n_pooled_qc_pos')} pooled QCs (POS)")
    roles = list(ROLE)
    axes[1].legend([Rectangle((0, 0), 1, 1, facecolor=ROLE[k], edgecolor="#333333", lw=0.3)
                    for k in roles], [ROLE_LABEL[k] for k in roles], ncol=6,
                   loc="upper center", bbox_to_anchor=(0.5, -0.9), handlelength=1.0,
                   columnspacing=1.2)


def qc_rsd(fig, gs, tag, compact=False):
    """Peak-picked only (slate) vs including gap-filled values (grey)."""
    set_panel("qc_rsd")
    vals, summ = D("qc_rsd_values"), D("qc_rsd_summary")
    axes = [fig.add_subplot(gs[i]) for i in range(2)]
    for ax, pol in zip(axes, POLS):
        thr = N(f"qc_rsd_threshold_pct_{pol}")
        for mode, colour, lab in [("peak_picked_only", SOIL, "detected peaks only"),
                                  ("including_gap_filled", "#A6A6A6",
                                   "incl. gap-filled values")]:
            v = vals[(vals.pol == pol) & (vals["mode"] == mode)].rsd_pct.to_numpy()
            s = summ[(summ.pol == pol) & (summ["mode"] == mode)].iloc[0]
            n = int(s.n_eligible)
            ax.step(v, np.arange(1, n + 1) / n, where="post", color=colour, lw=1.2,
                    label=lab if compact else f"{lab}\n{n:,} features, {s.pct_le_30} % ≤ {thr:g} %")
        ax.axvline(thr, color="black", lw=0.5, ls=":")
        ax.set_xlim(0, 150)
        ax.set_ylim(0, 1.0)
        ax.set_xlabel("pooled-QC RSD (%)")
        ax.text(0.97, 0.97, POL_LABEL[pol], transform=ax.transAxes, ha="right", va="top",
                fontsize=6.5, fontweight="bold", color="#555555")
        ax.legend(loc="lower right", fontsize=5, handlelength=1.2, borderaxespad=0.1,
                  labelspacing=0.6)
    axes[1].tick_params(labelleft=False)
    axes[0].set_ylabel("cumulative fraction of features")
    finding(axes[0], tag, f"Gap filling adds {N('qc_gapfill_feature_ratio_pos'):.0f}× more "
            f"features; only {N('qc_rsd_pct_le_threshold_including_gap_filled_pos'):.0f} % "
            "reproducible (POS)", x=-0.02)


def qc_internal_standards(fig, gs, tag):
    set_panel("qc_is")
    pts, info, seq = D("qc_is_points"), D("qc_is_info").set_index("pol"), D("qc_sequence")
    axes = np.array([[fig.add_subplot(gs[r, c]) for c in range(2)] for r in range(2)])
    for c, pol in enumerate(POLS):
        ipi, i = pts[pts.pol == pol], info.loc[pol]
        rt_tol, ppm_tol = i.rt_tol_min, i.ppm_tol
        spans, _ = _spans(pol)
        for r, (col, tol, ylab, lab) in enumerate([
                ("rt_minus_block_median", rt_tol, "RT − block median (min)",
                 f"±{rt_tol:g} min; block 2 − block 1: {i.block_shift_min:+.3f} min (median, "
                 f"n = {int(i.n_is_shift)} IS)"),
                ("ppm_error", ppm_tol, "mass error (ppm)", f"±{ppm_tol:g} ppm")]):
            ax = axes[r, c]
            _shade(ax, spans)
            ax.scatter(ipi.injection_order, ipi[col], s=4, color=POL[pol], alpha=0.7, lw=0,
                       zorder=3)
            for s_ in (-1, 1):
                ax.axhline(s_ * tol, color="black", lw=0.6, ls="--", zorder=2)
            ax.text(0.99, 0.97, lab, transform=ax.transAxes, ha="right", va="top", fontsize=5.8,
                    bbox=dict(facecolor="white", edgecolor="none", pad=0.5, alpha=0.8), zorder=4)
            ax.set_xlim(0.4, seq[seq.pol == pol].injection_order.max() + 0.6)
            lim = max(np.nanmax(np.abs(ipi[col])), tol) * 1.35
            ax.set_ylim(-lim, lim)
            if c == 0:
                ax.set_ylabel(ylab)
            if r == 0:
                ax.set_title(f"{POL_LABEL[pol]} ({int(i.n_is_detected)} IS detected)",
                             loc="left", fontsize=7, color=POL[pol], pad=2)
                ax.tick_params(labelbottom=False)
            else:
                ax.set_xlabel("injection order")
    letter(axes[0, 0], tag, x=-0.13, y=1.12,
           title=f"Internal standards: {N('qc_is_pct_rt_within'):.0f} % of RT values within "
                 f"±{N('qc_rt_tolerance_min_neg'):g} min of the block median, "
                 f"{N('qc_is_pct_ppm_within'):.0f} % of mass errors within "
                 f"±{N('qc_mass_tolerance_ppm_neg'):g} ppm (all injections)")


def qc_is_vs_cal(fig, gs, tag):
    set_panel("qc_is_vs_cal")
    ax = fig.add_subplot(gs)
    d = D("qc_is_vs_cal")
    for i, r in d.iterrows():
        ax.plot([r.cal_min_ng_per_ul, r.cal_max_ng_per_ul], [i, i], color="#AAAAAA", lw=3,
                solid_capstyle="butt")
        ax.scatter(r.soil_median_ng_per_ul, i, color=POL[r.pol], s=12, zorder=3)
    ax.set_yticks(range(len(d)))
    ax.set_yticklabels([f"{c.split('_')[0]} ({POL_LABEL[p]}, {n}/{m})" for c, p, n, m in
                        zip(d.compound, d.pol, d.n_samples_detected, d.n_samples)])
    ax.set_xscale("log")
    ax.set_xlabel(r"ng µL$^{-1}$ in extract (grey = calibration range)")
    ax.invert_yaxis()
    k = N("qc_is_n_within_cal")
    finding(ax, tag, f"{k} of {N('qc_is_n_detected')} detected standards "
            f"{'lies' if k == 1 else 'lie'} within the calibration range", x=-0.02, y=1.02)


def qc_recovery(ax, tag):
    set_panel("qc_recovery")
    d = D("qc_recovery")
    for i, r in d.iterrows():
        ax.plot([i, i], [r.R_app_min_pct, r.R_app_max_pct], color=POL[r.pol], lw=0.8)
        ax.scatter(i, r.R_app_median_pct, color=POL[r.pol], s=12, zorder=3)
    ax.axhline(100, color="black", lw=0.6, ls="--")
    ax.set_xticks(range(len(d)))
    ax.set_xticklabels([f"{c.split('_')[0]} {POL_LABEL[p]}" for c, p in zip(d.compound, d.pol)],
                       rotation=90)
    ax.set_ylabel("apparent recovery (%), median and range")
    ax.set_ylim(0, None)
    k, n = N("qc_recovery_n_above_100"), N("qc_recovery_n_standards")
    finding(ax, tag, f"Above 100 % for {'all ' + str(n) if k == n else f'{k} of {n}'} "
            "standards", x=-0.02, y=1.03)


# =========================================================================== annotation panels
def ann_upset(fig, gs, pol, first, tag):
    sets = D("ann_upset_sets")
    sets = sets[sets.pol == pol].sort_values("order")
    inter = D("ann_upset_intersections")
    inter = inter[inter.pol == pol].reset_index(drop=True)
    names = list(sets.source)
    sub = gs.subgridspec(2, 3, width_ratios=[0.22, 0.42, 1], height_ratios=[1, 0.62],
                         wspace=0.02, hspace=0.06)
    ax_bar = fig.add_subplot(sub[0, 2])
    ax_mat = fig.add_subplot(sub[1, 2], sharex=ax_bar)
    ax_set = fig.add_subplot(sub[1, 0], sharey=ax_mat)
    ax_lab = fig.add_subplot(sub[1, 1], sharey=ax_mat)
    x = np.arange(len(inter))
    ax_bar.bar(x, inter.n_features, color=POL[pol], width=0.7)
    for xi, n in zip(x, inter.n_features):
        ax_bar.text(xi, n, f"{n:,}", ha="center", va="bottom", fontsize=5, rotation=90)
    ax_bar.set_ylabel("features")
    ax_bar.tick_params(labelbottom=False, bottom=False)
    ax_bar.spines["bottom"].set_visible(False)
    ax_bar.set_ylim(0, inter.n_features.max() * 1.3)
    ax_bar.set_title(f"{POL_LABEL[pol]}: {N(f'ann_features_with_source_{pol}'):,} features "
                     "with a source", loc="left", fontsize=7, color=POL[pol], pad=2)
    ypos = {n: i for i, n in enumerate(names)}
    for xi, combo in zip(x, inter.intersection):
        ys = [ypos[m.strip()] for m in combo.split(";")]
        ax_mat.scatter([xi] * len(names), range(len(names)), s=10, color="#DDDDDD", lw=0)
        ax_mat.scatter([xi] * len(ys), ys, s=10, color="#222222", lw=0, zorder=3)
        if len(ys) > 1:
            ax_mat.plot([xi, xi], [min(ys), max(ys)], color="#222222", lw=0.8)
    ax_mat.set_axis_off()
    ax_mat.set_ylim(len(names) - 0.5, -0.5)
    sizes = sets.set_index("source").n_features.reindex(names)
    ax_set.barh(range(len(names)), sizes, color="#777777", height=0.6)
    for i, n in enumerate(sizes):
        ax_set.text(n, i, f"{n:,} ", ha="right", va="center", fontsize=5.5)
    ax_set.set_xlim(sizes.max() * 1.7, 0)
    ax_set.set_axis_off()
    ax_lab.set_axis_off()
    for n, i in ypos.items():
        ax_lab.text(0.04, i, n, va="center", ha="left", fontsize=6,
                    transform=ax_lab.get_yaxis_transform())
    if first:
        letter(ax_bar, tag, x=-0.62, y=1.2,
               title=f"{N('ann_single_source_pct_pos'):.0f} % of annotated POS features are "
                     "named by a single source")


def ann_confidence(ax, tag):
    set_panel("ann_confidence")
    conf = D("ann_confidence")
    ys = []
    for i, pol in enumerate(POLS):
        c = conf[conf.pol == pol]
        left = 0
        for _, r in c.iterrows():
            n = int(r.n_features)
            ax.barh(i, n, left=left, color=CONF[r.level], height=0.6, edgecolor="white", lw=0.4)
            ax.text(left + n / 2, i, f"{n:,}", ha="center", va="center", fontsize=5.5,
                    color="black" if r.level != "Bronze" else "white")
            left += n
        ax.text(left, i, f"  {left:,} of {int(c.n_all_features.iloc[0]):,}", va="center",
                fontsize=5.5)
        ys.append(i)
    ax.set_yticks(ys)
    ax.set_yticklabels([POL_LABEL[p] for p in POLS])
    for t, p in zip(ax.get_yticklabels(), POLS):
        t.set_color(POL[p])
    ax.invert_yaxis()
    ax.set_xlabel("annotated features")
    ax.set_xlim(0, ax.get_xlim()[1] * 1.3)
    handles = [plt.Rectangle((0, 0), 1, 1, color=CONF[k]) for k in ("Gold", "Silver", "Bronze")]
    ax.legend(handles, ["Gold", "Silver", "Bronze"], ncol=3, loc="lower right",
              bbox_to_anchor=(1.0, 1.0), handlelength=1, columnspacing=0.8)
    letter(ax, tag, x=-0.15, y=1.15,
           title=f"Most annotations are {N('ann_top_level_pos')} "
                 f"({N('ann_top_level_share_pos'):.0f} % of annotated POS features)")


def ann_classes(fig, gs, tag):
    set_panel("ann_classes")
    top = D("ann_classes_top")
    axes = [fig.add_subplot(gs[i]) for i in range(2)]
    for ax, pol in zip(axes, POLS):
        c = top[top.pol == pol]
        x = np.arange(len(c))
        ax.bar(x - 0.2, c.n_features, width=0.4, color=POL[pol], label="features")
        ax.bar(x + 0.2, c.n_unique_lipids, width=0.4, color="white", edgecolor=POL[pol],
               lw=0.7, label="unique lipids")
        ax.set_xticks(x)
        ax.set_xticklabels(c.lipid_class, rotation=90)
        ax.set_title(f"{POL_LABEL[pol]}: {len(c)} largest classes", loc="left", fontsize=7,
                     color=POL[pol], pad=2)
        ax.legend(loc="upper right", handlelength=1)
    axes[0].set_ylabel("Gold + Silver + Bronze")
    letter(axes[0], tag, x=-0.12, y=1.1,
           title=f"{N('ann_top3_classes_pos')} lead in positive mode, {N('ann_top_class_neg')} "
                 f"in negative mode; collapsing adducts removes "
                 f"{N('ann_adduct_redundancy_pct_pos'):.0f} % of POS features")


def ann_rt(fig, gs, tag):
    set_panel("ann_rt")
    ax = fig.add_subplot(gs)
    d = D("ann_rt_pc")
    lo, hi = np.inf, -np.inf
    notes = []
    for pol in POLS:
        t = d[d.pol == pol]
        ok, out = t[t.rt_check == "ok"], t[t.rt_check == "outlier"]
        ax.scatter(ok.rt_expected, ok.rt, s=5, color=POL[pol], lw=0, alpha=0.8,
                   label=f"{POL_LABEL[pol]} within tolerance ({len(ok)})")
        ax.scatter(out.rt_expected, out.rt, s=9, marker="x", color=POL[pol], lw=0.6,
                   label=f"{POL_LABEL[pol]} outlier ({len(out)})")
        lo, hi = min(lo, t.rt.min(), t.rt_expected.min()), max(hi, t.rt.max(), t.rt_expected.max())
        notes.append(f"{POL_LABEL[pol]} {N(f'ann_rt_untested_{pol}')}")
    ax.plot([lo, hi], [lo, hi], color="black", lw=0.6, ls="--")
    ax.set_xlabel("fitted RT (min)")
    ax.set_ylabel("observed RT (min)")
    ax.legend(loc="upper left", fontsize=5.5, handletextpad=0.2, borderaxespad=0.1)
    ax.text(0.02, 0.5, "untested, not drawn:\n" + ", ".join(notes), transform=ax.transAxes,
            ha="left", fontsize=5.5)
    letter(ax, tag, x=-0.15, y=1.04,
           title=f"RT model flags {N('ann_rt_outliers_pos')} of {N('ann_rt_tested_pos')} "
                 "tested PC names (POS)")


# =========================================================================== response panels
ADDUCT = {  # adduct -> (polarity, marker, label)
    "M+H": ("pos", "o", r"[M+H]$^{+}$"), "M+NH4": ("pos", "s", r"[M+NH$_4$]$^{+}$"),
    "M-H": ("neg", "o", r"[M−H]$^{-}$"), "M+HCOO": ("neg", "s", r"[M+HCOO]$^{-}$"),
}


def resp_rie(ax, tag):
    set_panel("resp_rie")
    d = D("resp_rie")
    order = d.drop_duplicates("class").sort_values("class_rank")["class"].tolist()
    rng = np.random.default_rng(0)  # horizontal jitter only, for legibility
    for add, (pol, mk, lab) in ADDUCT.items():
        s = d[d.adduct == add]
        x = s.class_rank + rng.uniform(-0.18, 0.18, len(s))
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
    letter(ax, tag, x=-0.1, y=1.03,
           title=f"Response spans {N('resp_rie_orders'):.0f} orders of magnitude across classes "
                 f"and adducts ({N('resp_n_standards')} standards)")


def resp_ie(ax, tag):
    set_panel("resp_ie")
    p = D("resp_ie_loco")
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
            f"{N('resp_ie_model').replace('_', ' ')}\nR² = {N('resp_ie_r2'):.2f}\n"
            f"RMSE = {N('resp_ie_rmse_log10'):.2f} log$_{{10}}$ "
            f"(~{N('resp_ie_fold_error'):.0f}-fold)\n"
            f"{N('resp_ie_n_standards')} POS standards, {N('resp_ie_n_classes')} classes",
            transform=ax.transAxes, va="top", fontsize=6)
    letter(ax, tag, x=-0.2, y=1.03,
           title=f"Response is not predictable across classes (R² = {N('resp_ie_r2'):.2f})")


def resp_nist(ax, tag, include_q3):
    set_panel("resp_nist")
    d = D("resp_nist")
    lo, hi = np.inf, -np.inf
    for pol in POLS:
        q1 = d[(d.pol == pol) & (d.route_group == "Q1")]
        lab = f"{POL_LABEL[pol]} Q1, n = {len(q1)}"
        if len(q1):
            lab += (f"\nmedian ratio {N(f'nist_q1_median_ratio_{pol}'):.2f}, "
                    f"{N(f'nist_q1_pct_within_2fold_{pol}'):.0f}% within 2-fold")
        ax.scatter(q1.medm, q1.estimate_nmol_per_ml, s=10, color=POL[pol], lw=0, zorder=3,
                   label=lab)
        sets = [q1]
        if include_q3:
            q3 = d[(d.pol == pol) & (d.route_group == "Q3")]
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
    letter(ax, tag, x=-0.2, y=1.03,
           title=f"Q1 vs NIST: median ratio {N('nist_q1_median_ratio_pos'):.2f}, "
                 f"{N('nist_q1_pct_within_2fold_pos'):.0f} % within 2-fold (POS)")


def route_tree(ax, tag):
    """Q1-Q4 decision tree; labels follow lib/soillifeatlas/quantify.py (option A).
    Schematic, no data; Q2 drawn but not used in this run."""
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
    ax.plot([xs[-1] + dw, 99.5, 99.5, xs[0]], [yq, yq, 2.5, 2.5], color="#AAAAAA", lw=0.6,
            ls=(0, (2, 1.5)))
    ax.annotate("", xy=(xs[0], yr - 6.2), xytext=(xs[0], 2.5),
                arrowprops=dict(arrowstyle="-|>", lw=0.6, color="#AAAAAA", mutation_scale=6))
    ax.text(99, yq + 1.2, "no", ha="right", va="bottom", fontsize=5, color="#555555")
    letter(ax, tag, x=0.0, y=1.0, title="Each feature gets one route; Q1–Q3 give an amount, Q4 "
           "keeps relative signal")


def quant_class_amounts(ax, tag, top=10):
    set_panel("quant_class_amounts")
    s = D("quant_class_amounts")
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


# =========================================================================== workflow panels
REPORT_HTML = "09_report/report_pos.html"
REPORT_PNG = "figures/rendered/report_pos.png"
REPORT_VIEW_CSS = (560, 300)  # narrow view: title and the summary tiles only
CHROME = "/Applications/Google Chrome.app/Contents/MacOS/Google Chrome"
# Stage, tool and Atlas-layer text, taken over from the 2026-11-24 draft schematic
# (figures/schematics/source/Fig2_workflow_2026-11-24.png). Schematic labels, not data.
STAGES = [
    # (number, title, what it does, tools, output, in Nextflow?)
    (1, "Convert", "raw → indexed mzML", "ThermoRawFileParser", None, True),
    (2, "Detect", "features per polarity;\npeak-picked and\ngap-filled kept apart", "MZmine 4.9",
     None, True),
    (3, "Check quality", "pooled QCs, internal\nstandards, blanks,\ncalibration, NIST",
     "workflow code", "QC report", True),
    (4, "Annotate", "lipid names with\nconfidence levels", "LipidSearch, MS-DIAL,\nGNPS, SIRIUS, "
     "MZmine", "annotated features", True),
    (5, "Amounts", "routes Q1–Q4,\nper g dry soil", "internal standards", "amounts per g soil",
     False),
    (6, "Report", "one-page report\nand result tables", "workflow code", "report + tables",
     False),
    (7, "Interpret", "spectral match to\nreference organisms", "Atlas matching", "Atlas matches",
     False),
]
ATLAS_LINKS = [(4, "identity\nreference spectra and names", False),
               (5, "response\nrelative ionisation by class", True),
               (7, "identity + occurrence\nwhere each lipid is found", True)]


def render_report() -> None:
    """Re-render the report screenshot when Chrome is available; otherwise keep the cached one."""
    html = input_path(REPORT_HTML)
    out = FIG_DIR / "rendered" / "report_pos.png"
    out.parent.mkdir(exist_ok=True)
    chrome = CHROME if os.path.exists(CHROME) else shutil.which("google-chrome")
    if chrome:
        w, h = REPORT_VIEW_CSS
        subprocess.run([chrome, "--headless=new", "--disable-gpu", "--hide-scrollbars",
                        "--force-device-scale-factor=2", f"--window-size={w},{h}",
                        f"--screenshot={out}", html.as_uri()],
                       check=True, capture_output=True)
    elif not out.exists():
        raise FileNotFoundError("no Chrome to render the report and no cached screenshot")


def report(ax):
    set_panel("report")
    render_report()
    show_image(ax, input_path(REPORT_PNG))
    ax.set_axis_on()
    ax.set_xticks([])
    ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(True)
        s.set_color("#999999")
        s.set_linewidth(0.5)


def workflow(ax):
    """Seven stages left to right; the Atlas as a band underneath feeding stages 4, 5 and 7."""
    set_panel("schematic_native")
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 40)
    ax.set_axis_off()
    x0, w, gap = 10.5, 11.3, 1.4
    top, h = 39, 22
    centres = {}
    ax.add_patch(FancyBboxPatch((0.3, top - h), 8.6, h, boxstyle="round,pad=0,rounding_size=1.2",
                                facecolor="#F2F2F2", edgecolor="#BBBBBB", lw=0.6))
    ax.text(4.6, top - 2.5, "Inputs", ha="center", va="top", fontsize=6.3, fontweight="bold")
    ax.text(4.6, top - 7, "raw files\n\nsample\nsheet\n\nAtlas\nrelease", ha="center", va="top",
            fontsize=5.5, linespacing=1.1)
    for i, (n, title, what, tools, out, nf) in enumerate(STAGES):
        x = x0 + i * (w + gap)
        centres[n] = x + w / 2
        ax.add_patch(FancyBboxPatch((x, top - h), w, h, boxstyle="round,pad=0,rounding_size=1.2",
                                    facecolor="white", edgecolor=SOIL, lw=0.9,
                                    linestyle="-" if nf else (0, (3, 2))))
        ax.add_patch(FancyBboxPatch((x, top - 5.2), w, 5.2,
                                    boxstyle="round,pad=0,rounding_size=1.2",
                                    facecolor=SOIL if nf else "#8A96A3", edgecolor="none"))
        ax.text(x + w / 2, top - 2.6, f"{n}  {title}", ha="center", va="center", fontsize=6,
                fontweight="bold", color="white")
        ax.text(x + w / 2, top - 7.2, what, ha="center", va="top", fontsize=5.3,
                linespacing=1.15)
        if i:
            ax.annotate("", xy=(x - 0.1, top - h / 2), xytext=(x - gap + 0.1, top - h / 2),
                        arrowprops=dict(arrowstyle="-|>", lw=0.6, color="#777777",
                                        mutation_scale=6))
        if out:
            ax.text(x + w / 2, top - h + 2.2, f"→ {out}", ha="center", va="bottom",
                    fontsize=5.2, color=SOIL, fontweight="bold")
    ax.annotate("", xy=(x0 - 0.1, top - h / 2), xytext=(9.0, top - h / 2),
                arrowprops=dict(arrowstyle="-|>", lw=0.6, color="#777777", mutation_scale=6))
    band_y, band_h = 0.8, 7.5
    ax.add_patch(FancyBboxPatch((x0, band_y), 7 * w + 6 * gap, band_h,
                                boxstyle="round,pad=0,rounding_size=1.2", facecolor=ATLAS,
                                alpha=0.12, edgecolor=ATLAS, lw=0.8))
    ax.text(x0 - 1.0, band_y + band_h / 2, "Soil Life\nAtlas", ha="right", va="center",
            fontsize=6.5, fontweight="bold", color=ATLAS)
    for stage, text, default in ATLAS_LINKS:
        cx = centres[stage]
        ax.text(cx, band_y + band_h / 2, text, ha="center", va="center", fontsize=5,
                color="#0B5F57", linespacing=1.15,
                bbox=dict(boxstyle="round,pad=0.3", facecolor="white", edgecolor=ATLAS, lw=0.6))
        ax.add_patch(FancyArrowPatch((cx, band_y + band_h + 0.2), (cx, top - h - 0.3),
                                     arrowstyle="-|>", mutation_scale=7, lw=0.9, color=ATLAS,
                                     linestyle="-" if default else (0, (2, 1.5))))


# =========================================================================== Atlas panels
def atlas_coverage(ax, tag):
    set_panel("atlas_coverage")
    cov = D("atlas_coverage")
    n_ms2 = N("atlas_n_ms2")
    colours = {"matched, no name": (ATLAS, "white"), "matched and named": (NAMED, "white"),
               "no Atlas match": (NOMATCH, "#333333")}
    left = 0
    for _, r in cov.iterrows():
        n = int(r.n)
        c, tc = colours[r.segment]
        ax.barh(0, n, left=left, color=c, height=0.6, edgecolor="white", lw=0.8)
        ax.text(left + n / 2, 0, f"{n:,}", ha="center", va="center", color=tc, fontsize=7,
                fontweight="bold")
        ax.text(left + n / 2, 0.45, f"{r.segment}\n{n / n_ms2 * 100:.0f} %", ha="center",
                va="bottom", fontsize=6)
        left += n
    ax.set_xlim(0, n_ms2)
    ax.set_ylim(-0.45, 1.25)
    ax.set_axis_off()
    finding(ax, tag, f"{N('atlas_pct_matched'):.0f} % of soil spectra match an Atlas reference; "
            f"{N('atlas_pct_unnamed_of_matched'):.0f} % of these have no name", x=0.0, y=1.08)


def atlas_breadth(ax, tag):
    set_panel("atlas_breadth")
    b = D("atlas_breadth")
    ax.bar(b.n_groups, b.unnamed, color=ATLAS, width=0.75, label="no name")
    ax.bar(b.n_groups, b.named, bottom=b.unnamed, color=NAMED, width=0.75, label="named")
    ax.set_xticks(b.n_groups)
    ax.set_xlabel("organism groups in which the lipid occurs")
    ax.set_ylabel("matched soil features")
    ax.legend(loc="upper right", handlelength=1)
    finding(ax, tag, f"{N('atlas_pct_single_group'):.0f} % point to a single organism group",
            x=-0.02)


def atlas_per_group(ax, tag):
    set_panel("atlas_per_group")
    g = D("atlas_per_group").set_index("group")
    order = g.occurs_pct.sort_values().index
    y = np.arange(len(order))
    for yi, grp in zip(y, order):
        ax.barh(yi, g.occurs_pct[grp], color=GROUP_COLOR[grp], alpha=0.35, height=0.7)
        ax.barh(yi, g.only_pct[grp], color=GROUP_COLOR[grp], height=0.7)
    ax.set_yticks(y)
    ax.set_yticklabels([f"{grp} ({int(g.ref_samples[grp])})" for grp in order])
    ax.tick_params(axis="y", length=0)
    ax.set_xlabel("matched soil features (%)")
    ax.text(0.98, 0.04, "pale: occurs in group\nsolid: only in this group\n(n reference samples)",
            transform=ax.transAxes, ha="right", va="bottom", fontsize=5.5, color="#444444")
    finding(ax, tag, "Matches point to all six organism groups", x=-0.02)


def atlas_class_group(ax, tag):
    set_panel("atlas_class_group")
    cg = D("atlas_class_group")
    order = cg.drop_duplicates("lipid_class").sort_values("class_rank")
    for _, r in cg.iterrows():
        ax.scatter(GROUP_ORDER.index(r.group), int(r.class_rank), s=r.share_pct * 1.1,
                   color=GROUP_COLOR[r.group], lw=0)
    ax.set_xticks(range(len(GROUP_ORDER)))
    ax.set_xticklabels(GROUP_ORDER, rotation=35, ha="right")
    ax.set_yticks(range(len(order)))
    ax.set_yticklabels([f"{r.lipid_class} ({int(r.n_features)})" for _, r in order.iterrows()])
    ax.set_xlim(-0.6, len(GROUP_ORDER) - 0.4)
    ax.set_ylim(len(order) - 0.5, -0.5)
    ax.tick_params(length=0)
    for s in ax.spines.values():
        s.set_visible(False)
    for v in (25, 50, 100):
        ax.scatter([], [], s=v * 1.1, color="#999999", lw=0, label=f"{v} %")
    ax.legend(title="named matches\nin group", loc="upper left", bbox_to_anchor=(1.0, 1.0),
              labelspacing=1.1, fontsize=5.5, title_fontsize=5.5, borderaxespad=0)
    finding(ax, tag, "Lipid classes differ in where they occur", x=-0.02)


def match_settings() -> str:
    return (f"modified cosine ≥ {N('atlas_match_min_cosine')}, ≥ {N('atlas_match_min_peaks')} "
            f"shared fragments, precursor ± {N('atlas_match_ppm'):g} ppm")


# =========================================================================== interpretation panels
def mirror(ax, role, tag, headline):
    """Soil MS/MS spectrum (up) against its Atlas reference (down), zoomed to the precursor
    range; the strongest shared fragments are labelled with their m/z."""
    set_panel(f"mirror_{role}")
    exs = D("interp_examples")
    ex = exs[exs.role == role].iloc[0]
    m = D(f"interp_mirror_{role}")
    soil, atl = m[m.side == "soil"], m[m.side == "atlas"]
    xmax = N(f"interp_xmax_{role}")
    ax.vlines(soil.mz, 0, soil.rel_intensity, color=SOIL, lw=0.7)  # rel_intensity is % of base
    ax.vlines(atl.mz, 0, -atl.rel_intensity, color=ATLAS, lw=0.7)
    ax.axhline(0, color="black", lw=0.5)
    for _, r in D(f"interp_labels_{role}").iterrows():
        ax.text(r.mz, r.height + 4, f"{r.mz:.2f}", ha="center", va="bottom", fontsize=5.5)
    ax.set_xlim(40, xmax)
    ax.set_ylim(-118, 145)
    ax.set_yticks([-100, 0, 100])
    ax.set_yticklabels(["100", "0", "100"])
    ax.set_xlabel("m/z")
    ax.set_ylabel("rel. intensity (%)")
    ax.text(0.99, 0.98, f"soil  ·  m/z {ex.soil_mz:.3f}, RT {ex.rt:.1f} min",
            transform=ax.transAxes, va="top", ha="right", fontsize=6, color=SOIL,
            fontweight="bold")
    ax.text(0.99, 0.03, f"Atlas reference  ·  m/z {ex.ref_mz:.3f}", transform=ax.transAxes,
            va="bottom", ha="right", fontsize=6, color=ATLAS, fontweight="bold")
    note = f"cosine {ex.cosine:.3f}\n{int(ex.matched_peaks)} shared fragments"
    ax.text(0.99, 0.36, note, transform=ax.transAxes, va="top", ha="right", fontsize=5.8,
            color="#333333")
    finding(ax, tag, headline.format(name=ex.lipid_name, phyla=_phyla(int(ex.n_ref_phyla))),
            x=-0.02)


def occurrence_matrix(ax, tag):
    """Both example lipids across every analysed phylum, phyla grouped by organism group."""
    set_panel("phylum_detection")
    occ, exs = D("interp_occurrence"), D("interp_examples").set_index("role")
    roles = ["broad", "restricted"]
    cols = occ[occ.role == "broad"].sort_values("col")
    for i, role in enumerate(roles):
        for _, r in occ[occ.role == role].iterrows():
            if r.detected:
                ax.scatter(r.col, i, s=20 + r.detection_rate * 130, color=GROUP_COLOR[r.group],
                           lw=0, zorder=3)
            else:
                ax.scatter(r.col, i, s=9, facecolor="white", edgecolor="#BBBBBB", lw=0.5,
                           zorder=3)
    start = 0
    for g in GROUP_ORDER:
        n = int((cols.group == g).sum())
        if not n:
            continue
        ax.add_patch(Rectangle((start - 0.45, -0.95), n - 0.1, 0.22, facecolor=GROUP_COLOR[g],
                               alpha=0.25, edgecolor="none"))
        ax.text(start + (n - 1) / 2, -1.02, g, ha="center", va="bottom", fontsize=6,
                fontweight="bold", color=GROUP_COLOR[g])
        start += n
    ax.set_xticks(range(len(cols)))
    ax.set_xticklabels(cols.phylum, rotation=45, ha="right", fontsize=5.8)
    ax.set_yticks(range(len(roles)))
    ax.set_yticklabels([f"{exs.loc[r, 'lipid_name']}\n({_phyla(int(exs.loc[r, 'n_ref_phyla']))})"
                        for r in roles])
    ax.set_xlim(-0.6, len(cols) - 0.4)
    ax.set_ylim(len(roles) - 0.35, -1.25)
    ax.tick_params(length=0)
    for s in ax.spines.values():
        s.set_visible(False)
    for v in (0.25, 0.5, 1.0):
        ax.scatter([], [], s=20 + v * 130, color="#999999", lw=0, label=f"{v * 100:.0f} %")
    ax.scatter([], [], s=9, facecolor="white", edgecolor="#BBBBBB", lw=0.5, label="not found")
    ax.legend(title="reference samples\nwith the lipid", loc="upper left",
              bbox_to_anchor=(1.0, 1.0), labelspacing=1.0, fontsize=5.5, title_fontsize=5.5,
              borderaxespad=0)
    parts = []
    for role in roles:
        gs = N(f"interp_groups_{role}").split(";")
        where = gs[0].lower() if len(gs) == 1 else f"{len(gs)} of {len(GROUP_ORDER)} groups"
        parts.append(f"{exs.loc[role, 'lipid_name']} occurs in {where}")
    finding(ax, tag, "; ".join(parts), x=-0.02, y=1.02)


# =========================================================================== application panels
TREAT = ["ambient", "future", "ambient_drought", "future_drought"]
TREAT_LABEL = {"ambient": "ambient", "future": "future", "ambient_drought": "ambient + drought",
               "future_drought": "future + drought"}
TREAT_COLOR = {"ambient": "#BBBBBB", "future": OI["orange"], "ambient_drought": OI["green"],
               "future_drought": OI["purple"]}


def app_design(ax, tag):
    set_panel("app_design")
    s = D("app_design")
    ax.set_xlim(0, 2)
    ax.set_ylim(0, 2)
    ax.set_axis_off()
    for i, wc in enumerate(["no", "yes"]):
        for j, dr in enumerate(["no", "yes"]):
            cell = s[(s.warming_co2 == wc) & (s.drought == dr)]
            t = cell.treatment.iloc[0]
            ax.add_patch(plt.Rectangle((i + 0.04, 1 - j + 0.04), 0.92, 0.92,
                                       facecolor=TREAT_COLOR[t], alpha=0.35, edgecolor="none"))
            lines = [f"{int((cell.used_pos == 'yes').sum())} plots"]
            ax.text(i + 0.5, 1 - j + 0.8, TREAT_LABEL[t].replace(" + ", "\n+ "), ha="center",
                    va="center", fontweight="bold", fontsize=6)
            ax.text(i + 0.5, 1 - j + 0.3, "\n".join(lines), ha="center", va="center",
                    fontsize=5.5)
    ax.text(0.5, 2.06, "ambient", ha="center", fontsize=6)
    ax.text(1.5, 2.06, "+3 °C, +300 ppm CO$_2$", ha="center", fontsize=6)
    ax.text(-0.04, 1.5, "no drought", ha="right", va="center", fontsize=6, rotation=90)
    ax.text(-0.04, 0.5, "drought", ha="right", va="center", fontsize=6, rotation=90)
    finding(ax, tag, f"{N('app_soils_pos')} ClimGrass soils, 2 × 2 design",
            x=-0.02, y=1.1)


def app_class_profiles(ax, tag, top=8):
    set_panel("app_class_profiles")
    p = D("app_class_profiles")
    order = p.groupby("lipid_class")["mean"].max().sort_values(ascending=False).head(top).index
    x = np.arange(len(order))
    w = 0.2
    for k, t in enumerate(TREAT):
        s = p[p.treatment == t].set_index("lipid_class").reindex(order)
        ax.bar(x + (k - 1.5) * w, s["mean"], w, yerr=s["std"], color=TREAT_COLOR[t],
               error_kw=dict(lw=0.4, capsize=1), label=TREAT_LABEL[t])
    ax.set_yscale("log")
    ax.set_xticks(x)
    ax.set_xticklabels(order, rotation=90)
    ax.set_ylabel(r"µg g$^{-1}$ (mean ± SD of plots)")
    ax.legend(fontsize=5.5, ncol=2)
    letter(ax, tag, x=-0.15, y=1.03, title="Class amounts by treatment (POS, common set)")


def app_responsive(ax, tag):
    set_panel("app_responsive")
    rows = D("app_responsive")
    y = np.arange(len(rows))
    ax.barh(y, rows.log2fc_drought_vs_ambient, color=[POL[p] for p in rows.pol], height=0.7)
    ax.axvline(0, color="black", lw=0.5)
    ax.set_yticks(y)
    ax.set_yticklabels(rows.lipid_name, fontsize=5)
    ax.set_xlabel(r"log$_2$ fold change, drought vs ambient")
    ax.legend([plt.Rectangle((0, 0), 1, 1, color=POL[p]) for p in POLS],
              [POL_LABEL[p] for p in POLS], loc="lower right", fontsize=5.5)
    finding(ax, tag, f"{N('app_n_responsive_pos')} POS and {N('app_n_responsive_neg')} NEG "
            "lipids change with drought", x=-0.45, y=1.02)


def app_occurrence(ax, tag):
    set_panel("app_occurrence")
    d = D("app_occurrence")
    y = np.arange(len(d))
    ax.barh(y, d.n_ref_phyla.fillna(0), color=ATLAS, height=0.7)
    for yy in y[d.n_ref_phyla.isna().to_numpy()]:
        ax.text(0.2, yy, "no Atlas match", va="center", fontsize=5, color="#777777")
    ax.set_yticks(y)
    ax.set_yticklabels(d.lipid_name, fontsize=5)
    ax.xaxis.set_major_locator(plt.MaxNLocator(integer=True))
    ax.set_xlabel("phyla in which the reference occurs")
    finding(ax, tag, f"{N('app_n_responsive_pos_matched')} of {N('app_n_responsive_pos')} POS "
            "lipids match an Atlas reference", x=-0.3, y=1.02)


def app_trace(fig, gs, tag, lipid):
    """One responsive lipid: amount per treatment (top) and where its reference occurs (bottom)."""
    set_panel("app_trace")
    t = D("app_trace")
    r = t[t.lipid_name == lipid].iloc[0]
    sub = gs.subgridspec(2, 1, height_ratios=[1, 0.62], hspace=0.9)
    ax = fig.add_subplot(sub[0])
    vals = [r[f"mean_ng_g_{x}"] for x in TREAT]
    ax.bar(range(4), vals, color=[TREAT_COLOR[x] for x in TREAT], width=0.7)
    ax.set_xticks(range(4), [TREAT_LABEL[x].replace(" + ", "\n+ ") for x in TREAT], fontsize=5.5)
    ax.set_ylabel(r"ng g$^{-1}$ dry soil")
    ax.set_ylim(0, max(vals) * 1.25)
    finding(ax, tag, f"{r.lipid_name}: lower under drought (q = {r.q_drought:.3f})", x=-0.02,
            y=1.06)
    g = D("app_trace_groups")
    g = g[g.lipid_name == lipid].set_index("group").reindex(GROUP_ORDER)
    bx = fig.add_subplot(sub[1])
    bx.set_axis_off()
    for k, grp in enumerate(GROUP_ORDER):
        n, tot = int(g.loc[grp, "n_phyla_with_ref"]), int(g.loc[grp, "n_phyla_analysed"])
        bx.text(-0.6, -k, grp, ha="right", va="center", fontsize=5.6)
        for m in range(tot):
            bx.scatter(m, -k, s=20, color=GROUP_COLOR[grp] if m < n else "white",
                       edgecolor=GROUP_COLOR[grp], lw=0.6, clip_on=False)
        bx.text(tot - 0.4, -k, f"{n} of {tot}", va="center", fontsize=5.2, color="#555555")
    bx.set_xlim(-4.2, 4.5)
    bx.set_ylim(-len(GROUP_ORDER) + 0.5, 0.5)
    bx.text(-4.2, 0.9, f"Atlas reference (cosine {r.cosine:.2f}) occurs in "
            f"{int(r.n_ref_phyla)} phyla", fontsize=6, va="bottom")


# =========================================================================== graphical abstract
def _box(ax, x, y, w, h, fc, ec="none", lw=0.8, ls="-", r=1.5, alpha=1.0):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle=f"round,pad=0,rounding_size={r}",
                                facecolor=fc, edgecolor=ec, lw=lw, linestyle=ls, alpha=alpha))


def graphical_abstract_panel(ax):
    set_panel("ga")
    ax.set_xlim(0, 200)
    ax.set_ylim(0, 80)
    ax.set_axis_off()
    ax.text(2, 77, "Which organisms could have made the lipids in a soil?", fontsize=10,
            fontweight="bold", va="top", color=SOIL)
    ax.text(2, 70.5, "The Soil Life Atlas links soil lipid spectra to reference organisms; an "
            "open workflow measures, checks and annotates the soil data", fontsize=6.8,
            va="top", color="#444444")
    # left: the Atlas
    per_group = D("atlas_release_groups").set_index("group").ref_samples
    _box(ax, 2, 4, 58, 58, ATLAS, alpha=0.10, ec=ATLAS, lw=0.9)
    ax.text(6, 58.5, "Soil Life Atlas", fontsize=8.5, fontweight="bold", color=ATLAS, va="top")
    ax.text(6, 53.5, f"{N('atlas_release_ref_samples')} reference samples in "
            f"{N('atlas_release_phyla')} analysed phyla\nspectra · response · occurrence", fontsize=5.6,
            color="#335", va="top", linespacing=1.3)
    ymax = per_group.max()
    for i, grp in enumerate(GROUP_ORDER):
        y = 41 - i * 6.2
        ax.scatter(8, y, s=40, color=GROUP_COLOR[grp], lw=0, zorder=3)
        ax.text(11, y, grp, fontsize=6.4, va="center")
        w = 22 * per_group[grp] / ymax
        _box(ax, 32, y - 1.6, w, 3.2, GROUP_COLOR[grp], r=0.8, alpha=0.85)
        ax.text(33 + w, y, f"{int(per_group[grp])}", fontsize=5.8, va="center")
    ax.add_patch(FancyArrowPatch((60.5, 20), (71.5, 13.4), arrowstyle="-|>", mutation_scale=9,
                                 lw=1.4, color=ATLAS, connectionstyle="arc3,rad=0.15"))
    ax.text(62, 22, "reference\nspectra", fontsize=5.2, color=ATLAS, va="bottom")
    # middle: the workflow
    ax.text(72, 58.5, "Soil sample", fontsize=8.5, fontweight="bold", color=SOIL, va="top")
    ax.text(72, 53.5, "LC–MS/MS, positive and negative mode", fontsize=5.6, color="#444",
            va="top")
    steps = [("detect features", "peak-picked, kept apart\nfrom gap-filled", True),
             ("check quality", "pooled QCs, standards,\nblanks", True),
             ("annotate", "names with confidence", True),
             ("match to the Atlas", "spectral similarity", False)]
    for i, (t, d, nf) in enumerate(steps):
        y = 44 - i * 10.2
        _box(ax, 72, y - 4, 44, 8, "white", ec=SOIL if nf else ATLAS, lw=0.9,
             ls="-" if nf else (0, (3, 2)))
        ax.text(75, y, t, fontsize=6.6, fontweight="bold", va="center",
                color=SOIL if nf else ATLAS)
        if i:
            ax.annotate("", xy=(94, y + 4.1), xytext=(94, y + 6.1),
                        arrowprops=dict(arrowstyle="-|>", lw=0.6, color="#888",
                                        mutation_scale=5))
    ax.text(72, 7.5, "open nf-core workflow · laptop or HPC", fontsize=5.4, color="#666")
    ax.add_patch(FancyArrowPatch((118, 33), (126.5, 33),
                                 arrowstyle="simple,head_width=4,head_length=4",
                                 mutation_scale=2.2, color=SOIL, lw=0))
    # right: what the profile gains
    n_ms2, n_match = N("atlas_n_ms2"), N("atlas_n_matched")
    ax.text(128, 58.5, "What the soil profile gains", fontsize=8.5, fontweight="bold",
            color=SOIL, va="top")
    ax.text(128, 53.5, "ClimGrass grassland soils, positive mode", fontsize=5.6, color="#444",
            va="top")
    x, W = 128, 68
    ax.text(x, 49.5, f"{N('atlas_pct_matched'):.0f} %", fontsize=14, fontweight="bold",
            color=ATLAS, va="top")
    ax.text(x + 18, 47, f"of soil MS/MS spectra match an Atlas\nreference ({n_match:,} of "
            f"{n_ms2:,})", fontsize=5.6, va="center", linespacing=1.2)
    colours = {"matched, no name": ATLAS, "matched and named": NAMED, "no Atlas match": NOMATCH}
    left = x
    for _, r in D("atlas_coverage").iterrows():
        w = W * r.n / n_ms2
        ax.add_patch(FancyBboxPatch((left, 37.5), w, 3.6, boxstyle="square,pad=0",
                                    facecolor=colours[r.segment], edgecolor="white", lw=0.6))
        left += w
    ax.text(x, 35.5, f"{N('atlas_pct_unnamed_of_matched'):.0f} % of the matched lipids have "
            "no name;\nthe Atlas still shows where they occur", fontsize=5.6, va="top",
            color="#0B5F57", linespacing=1.2)
    occ = D("atlas_per_group").set_index("group").occurs_pct
    ax.text(x, 27, "matched soil lipids found in each group", fontsize=5.8, va="top",
            color="#444")
    for i, grp in enumerate(GROUP_ORDER):
        y = 21.5 - i * 3.9
        pct = occ[grp]
        ax.text(x, y, grp, fontsize=5.6, va="center")
        _box(ax, x + 17, y - 1.2, 42 * pct / 100, 2.4, GROUP_COLOR[grp], r=0.6, alpha=0.9)
        ax.text(x + 18 + 42 * pct / 100, y, f"{pct:.0f} %", fontsize=5.4, va="center")


# =========================================================================== the Atlas (Fig 1)
LAYER_SHORT = {"consensus": "features measured", "quality_filtered": "quality-filtered",
               "enriched": "phylum-enriched lipids", "with_msms": "with a usable MS/MS spectrum",
               "annotated": "with an accepted name"}


def release_collection(ax, tag):
    """Reference samples per phylum, grouped and coloured by organism group."""
    set_panel("release")
    ph = D("release_phyla")
    ph["analysed"] = ph.analysed.astype(str) == "True"
    y, ticks, labels, gy = 0.0, [], [], []
    for grp in GROUP_ORDER:
        sub = ph[ph.group == grp]
        top = y
        for _, r in sub.iterrows():
            c = GROUP_COLOR[grp]
            ax.barh(y, r.samples_pos, color=c if r.analysed else "white", edgecolor=c,
                    lw=0.7, height=0.72, hatch=None if r.analysed else "////")
            ax.text(r.samples_pos + 0.6, y, f"{int(r.samples_pos)}", va="center", fontsize=5.8)
            ticks.append(y)
            labels.append(r.phylum)
            y += 1
        gy.append((grp, top, y - 1))
        y += 0.7
    ax.set_yticks(ticks, labels, fontsize=6)
    ax.set_ylim(y - 0.4, -0.8)
    ax.tick_params(axis="y", length=0)
    for grp, a, b in gy:
        ax.annotate("", xy=(-0.47, a - 0.35), xytext=(-0.47, b + 0.35),
                    xycoords=("axes fraction", "data"),
                    arrowprops=dict(arrowstyle="-", lw=2.2, color=GROUP_COLOR[grp]),
                    annotation_clip=False)
        ax.text(-0.5, (a + b) / 2, grp, transform=ax.get_yaxis_transform(), ha="right",
                va="center", fontsize=6.4, fontweight="bold", color="#333")
    ax.set_xlabel("reference samples (positive mode)")
    ax.spines["left"].set_visible(False)
    ax.text(0.98, 0.02, "hatched: collected, not yet analysed", transform=ax.transAxes,
            ha="right", va="bottom", fontsize=5.6, color="#555")
    finding(ax, tag, f"{N('release_reference_samples_pos')} reference samples from "
            f"{N('release_phyla_collected')} phyla in six organism groups", x=-0.52, y=1.02)


def release_layers(ax, tag):
    """From everything measured to the curated Atlas lipids, both polarities (log scale)."""
    set_panel("release")
    L = D("release_layers")
    ypos = np.arange(len(L))
    h = 0.36
    for off, pol, alpha in ((-h / 2, "pos", 1.0), (h / 2, "neg", 0.45)):
        v = L[pol].astype(float)
        ax.barh(ypos + off, v, height=h, color=ATLAS, alpha=alpha, lw=0,
                label=POL_LABEL.get(pol, pol.upper()))
        for yy, val in zip(ypos + off, v):
            ax.text(val * 1.15, yy, f"{int(val):,}", va="center", fontsize=5.6)
    ax.set_xscale("log")
    ax.set_xlim(1e3, 2e6)
    ax.set_yticks(ypos, [LAYER_SHORT[i] for i in L["id"]], fontsize=6.2)
    ax.set_ylim(len(L) - 0.4, -0.6)
    ax.tick_params(axis="y", length=0)
    ax.spines["left"].set_visible(False)
    ax.set_xlabel("features")
    enr = int(np.flatnonzero(L["id"].to_numpy() == "enriched")[0])
    ax.axhspan(enr - 0.48, enr + 0.48, color=ATLAS, alpha=0.08, lw=0)
    ax.legend(loc="lower right", handlelength=1, fontsize=5.8)
    finding(ax, tag, f"{N('release_lipids_pos'):,} POS and {N('release_lipids_neg'):,} NEG "
            "lipids make up the Atlas", x=-0.76, y=1.03)


def release_record(ax, tag):
    """What one Atlas record holds, and where the workflow uses it (schematic)."""
    set_panel("release_record")
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 40)
    ax.set_axis_off()
    items = [("Identity", "m/z, retention time,\nMS/MS spectrum, name", "annotation"),
             ("Occurrence", "phyla and organism groups\nin which it was detected",
              "interpretation"),
             ("Response", "class-level ionisation\nresponse", "amounts")]
    for i, (head, body, use) in enumerate(items):
        x = 1 + i * 33.3
        _box(ax, x, 13, 31, 23, ATLAS, alpha=0.10, ec=ATLAS, lw=0.8, r=1.2)
        ax.text(x + 2, 32.5, head, fontsize=7, fontweight="bold", color=ATLAS, va="top")
        ax.text(x + 2, 26, body, fontsize=5.8, va="top", linespacing=1.25)
        ax.text(x + 15.5, 7, f"used for {use}", fontsize=5.8, ha="center", va="center",
                color=SOIL, style="italic")
    ax.text(0, 42, tag, fontsize=9, fontweight="bold", va="bottom", ha="right")
    ax.text(1, 42, "Each Atlas lipid is a record that the workflow can use", fontsize=7,
            va="bottom")


# =========================================================================== figures
def graphical_abstract():
    P = Panels("GA", "graphical_abstract")
    fig = plt.figure(figsize=(DOUBLE, 76 * MM))
    graphical_abstract_panel(fig.add_axes([0.0, 0.0, 1.0, 1.0]))
    P.add("–", "REAL", "Atlas group counts from the release export (figures/inputs/); match "
          "numbers from 06_atlas_match/pos", ["ga"])
    return save(fig, "GA", "graphical_abstract"), P.rows


def fig1_atlas():
    P = Panels(1, "atlas")
    fig = plt.figure(figsize=(DOUBLE, 100 * MM))
    release_collection(fig.add_axes([0.2, 0.1, 0.25, 0.8]), "A")
    P.add("A", "REAL", f"{N('release_label')}; exported by export_groups.py", ["release"])
    release_layers(fig.add_axes([0.71, 0.47, 0.25, 0.43]), "B")
    P.add("B", "REAL", f"{N('release_label')}; feature layers of the release", ["release"])
    release_record(fig.add_axes([0.52, 0.03, 0.46, 0.3]), "C")
    P.add("C", "SCHEMATIC", "drawn natively; no data", ["release_record"])
    return save(fig, 1, "atlas"), P.rows


def fig2_workflow():
    P = Panels(2, "workflow")
    fig = plt.figure(figsize=(DOUBLE, 118 * MM))
    workflow(fig.add_axes([0.005, 0.53, 0.99, 0.43]))
    fig_letter(fig, "A", 0.005, 0.995, "Seven stages from raw files to Atlas-linked lipid "
               "profiles; the Atlas enters at annotation, amounts and interpretation")
    P.add("A", "SCHEMATIC", "drawn natively; stage and tool text from the 2026-11-24 draft",
          ["schematic_native"])
    qc_rsd(fig, fig.add_gridspec(1, 2, left=0.07, right=0.56, top=0.43, bottom=0.09,
                                 wspace=0.08), "B", compact=True)
    P.add("B", "REAL", "pooled-QC RSD per feature; eligibility rule in rsd_summary.tsv",
          ["qc_rsd"])
    report(fig.add_axes([0.62, 0.05, 0.36, 0.38]))
    fig_letter(fig, "C", 0.605, 0.485, "Every run ends in a one-page report")
    P.add("C", "REAL", "summary tiles of report_pos.html; the 'lipids with an amount' tile is "
          "stage-5 output", ["report"])
    return save(fig, 2, "workflow"), P.rows


def fig3_atlas():
    P = Panels(3, "atlas_matching")
    fig = plt.figure(figsize=(DOUBLE, 115 * MM))
    gs = fig.add_gridspec(2, 3, height_ratios=[0.3, 1], width_ratios=[0.9, 1, 1.1],
                          hspace=0.55, wspace=0.75, left=0.07, right=0.9, top=0.92, bottom=0.14)
    atlas_coverage(fig.add_subplot(gs[0, :]), "A")
    P.add("A", "REAL", "POS; " + match_settings(), ["atlas_coverage"])
    atlas_breadth(fig.add_subplot(gs[1, 0]), "B")
    P.add("B", "REAL", "organism groups from the release (figures/inputs/atlas_groups.csv)",
          ["atlas_breadth"])
    atlas_per_group(fig.add_subplot(gs[1, 1]), "C")
    P.add("C", "REAL", "% of matched features; reference samples per group in brackets "
          "(matching denominators)", ["atlas_per_group"])
    atlas_class_group(fig.add_subplot(gs[1, 2]), "D")
    P.add("D", "REAL", "8 most frequent classes among named matches", ["atlas_class_group"])
    return save(fig, 3, "atlas_matching"), P.rows


def fig4_interpretation():
    P = Panels(4, "interpretation")
    fig = plt.figure(figsize=(DOUBLE, 110 * MM))
    gs = fig.add_gridspec(2, 2, height_ratios=[1, 0.62], hspace=0.62, wspace=0.25, left=0.1,
                          right=0.86, top=0.92, bottom=0.18)
    mirror(fig.add_subplot(gs[0, 0]), "broad", "A",
           "{name}: soil spectrum matches a widespread reference ({phyla})")
    P.add("A", "REAL", "zoomed to m/z 40–precursor+25; example proposed, to confirm",
          ["mirror_broad"])
    mirror(fig.add_subplot(gs[0, 1]), "restricted", "B",
           "{name}: soil spectrum matches a restricted reference ({phyla})")
    P.add("B", "REAL", "zoomed; one reference peak above the range noted in the panel",
          ["mirror_restricted"])
    occurrence_matrix(fig.add_subplot(gs[1, :]), "C")
    P.add("C", "REAL", "all 16 analysed phyla, grouped by organism group", ["phylum_detection"])
    return save(fig, 4, "interpretation"), P.rows


def fig5_application():
    g = gates()
    panels = ["A", "B", "C"] + (["D"] if g.get("fig7E_lipid") else [])
    if not g["stage5_fix_confirmed"]:
        return _not_built(5, "application", "main", "waits for stage5_fix_confirmed", panels)
    P = Panels(5, "application")
    fig = plt.figure(figsize=(DOUBLE, 150 * MM))
    gs = fig.add_gridspec(2, 3, width_ratios=[1, 1.05, 1.05], height_ratios=[1, 1],
                          hspace=0.45, wspace=0.75, left=0.06, right=0.98, top=0.93,
                          bottom=0.1)
    app_design(fig.add_subplot(gs[0, 0]), "A")
    P.add("A", "REAL", "ClimGrass design; drying–rewetting validation only in the caption (#11)",
          ["app_design"])
    app_responsive(fig.add_subplot(gs[:, 1]), "B")
    P.add("B", "REAL", "drought-responsive lipids, q < 0.05", ["app_responsive"])
    app_occurrence(fig.add_subplot(gs[:, 2]), "C")
    P.add("C", "REAL", "Atlas occurrence of the POS responsive lipids", ["app_occurrence"])
    if len(panels) == 4:
        app_trace(fig, gs[1, 0], "D", g["fig7E_lipid"])
        P.add("D", "REAL", f"traced lipid: {g['fig7E_lipid']}", ["app_trace"])
    return save(fig, 5, "application"), P.rows


def figS1_quality():
    P = Panels("S1", "quality", "supplement")
    fig = plt.figure(figsize=(DOUBLE, 165 * MM))
    outer = fig.add_gridspec(3, 1, height_ratios=[0.55, 1.0, 1.35], hspace=0.62, left=0.09,
                             right=0.98, top=0.95, bottom=0.06)
    qc_sequence(fig, outer[0].subgridspec(2, 1, hspace=1.5), "A")
    P.add("A", "REAL", "POS and NEG sequences, both blocks", ["qc_sequence"])
    qc_rsd(fig, outer[1].subgridspec(1, 2, wspace=0.08), "B")
    P.add("B", "REAL", "ECDF, x cut at 150 %; n and % ≤ 30 % from rsd_summary.tsv", ["qc_rsd"])
    qc_internal_standards(fig, outer[2].subgridspec(2, 2, hspace=0.12, wspace=0.15), "C")
    P.add("C", "REAL", "RT relative to recorded block median; tolerances from provenance.yaml",
          ["qc_is"])
    return save(fig, "S1", "quality"), P.rows


def figS2_annotation():
    P = Panels("S2", "annotation", "supplement")
    fig = plt.figure(figsize=(DOUBLE, 185 * MM))
    outer = fig.add_gridspec(3, 1, height_ratios=[1.0, 1.0, 0.95], hspace=0.5, left=0.08,
                             right=0.98, top=0.93, bottom=0.07)
    set_panel("ann_upset")
    top = outer[0].subgridspec(1, 2, wspace=0.12)
    ann_upset(fig, top[0], "pos", True, "A")
    ann_upset(fig, top[1], "neg", False, "A")
    P.add("A", "REAL", "exclusive intersections, POS and NEG", ["ann_upset"])
    mid = outer[1].subgridspec(1, 2, width_ratios=[1.1, 1], wspace=0.3)
    ann_confidence(fig.add_subplot(mid[0]), "B")
    P.add("B", "REAL", "ClimGrass only", ["ann_confidence"])
    ann_rt(fig, mid[1], "D")
    P.add("D", "REAL", "PC; untested lipids counted, not drawn", ["ann_rt"])
    ann_classes(fig, outer[2].subgridspec(1, 2, wspace=0.15), "C")
    P.add("C", "REAL", "12 largest classes per polarity", ["ann_classes"])
    return save(fig, "S2", "annotation"), P.rows


def figS3_response():
    P = Panels("S3", "response_routes", "supplement")
    include_q3 = bool(gates()["stage5_fix_confirmed"])
    fig = plt.figure(figsize=(DOUBLE, 160 * MM))
    outer = fig.add_gridspec(2, 1, height_ratios=[1.0, 1.0], hspace=0.45, left=0.08,
                             right=0.98, top=0.95, bottom=0.08)
    r1 = outer[0].subgridspec(1, 2, width_ratios=[1.55, 1], wspace=0.28)
    resp_rie(fig.add_subplot(r1[0]), "A")
    P.add("A", "REAL", "RIE relative to LPE", ["resp_rie"])
    resp_ie(fig.add_subplot(r1[1]), "B")
    P.add("B", "REAL", "leave-one-class-out; statistics from ie_model_stats.tsv", ["resp_ie"])
    r2 = outer[1].subgridspec(1, 2, width_ratios=[1.45, 1], wspace=0.28)
    route_tree(fig.add_subplot(r2[0]), "C")
    P.add("C", "SCHEMATIC", "drawn natively; Q2 not used in this run", ["route_tree"])
    resp_nist(fig.add_subplot(r2[1]), "D", include_q3)
    P.add("D", "REAL", "Q1 and Q3" if include_q3 else "Q1 only; Q3 waits for "
          "stage5_fix_confirmed", ["resp_nist"])
    return save(fig, "S3", "response_routes"), P.rows


def figS4_recovery():
    if not gates()["reconstitution_volume_confirmed"]:
        return _not_built("S4", "recovery", "supplement",
                          "waits for reconstitution_volume_confirmed", ["A", "B"])
    P = Panels("S4", "recovery", "supplement")
    fig = plt.figure(figsize=(DOUBLE, 85 * MM))
    gs = fig.add_gridspec(1, 2, wspace=0.5, left=0.15, right=0.98, top=0.88, bottom=0.25)
    qc_is_vs_cal(fig, gs[0, 0], "A")
    P.add("A", "REAL", "IS in soil vs calibration range", ["qc_is_vs_cal"])
    qc_recovery(fig.add_subplot(gs[0, 1]), "B")
    P.add("B", "REAL", "apparent recovery", ["qc_recovery"])
    return save(fig, "S4", "recovery"), P.rows


def figS5_amounts():
    if not gates()["stage5_fix_confirmed"]:
        return _not_built("S5", "amounts", "supplement", "waits for stage5_fix_confirmed",
                          ["A", "B"])
    P = Panels("S5", "amounts", "supplement")
    fig = plt.figure(figsize=(DOUBLE, 85 * MM))
    gs = fig.add_gridspec(1, 2, wspace=0.3, left=0.08, right=0.98, top=0.88, bottom=0.2)
    quant_class_amounts(fig.add_subplot(gs[0, 0]), "A")
    P.add("A", "REAL", "class amounts, common lipid set", ["quant_class_amounts"])
    app_class_profiles(fig.add_subplot(gs[0, 1]), "B")
    P.add("B", "REAL", "class amounts by treatment (POS)", ["app_class_profiles"])
    return save(fig, "S5", "amounts"), P.rows


FIGURES = [graphical_abstract, fig1_atlas, fig2_workflow, fig3_atlas, fig4_interpretation, fig5_application,
           figS1_quality, figS2_annotation, figS3_response, figS4_recovery, figS5_amounts]


def main() -> int:
    for build in FIGURES:
        for p in build()[0]:
            print("wrote", p)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
