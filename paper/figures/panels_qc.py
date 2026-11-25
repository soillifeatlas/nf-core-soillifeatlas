"""Figure 3 - analytical quality (ClimGrass)."""
from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Rectangle

from common import (DOUBLE, MM, POL, POL_LABEL, ROLE, ROLE_LABEL, Panels, gates, letter, pending,
                    read_table, read_yaml, save, set_panel)

POLS = ("pos", "neg")


def block_spans(inj: pd.DataFrame) -> dict[int, tuple[int, int]]:
    g = inj.groupby("block").injection_order
    return {b: (int(v.min()), int(v.max())) for b, v in g}


def shade_blocks(ax, spans):
    for b, (lo, hi) in spans.items():
        if b % 2 == 1:
            ax.axvspan(lo - 0.5, hi + 0.5, color="#EFEFEF", zorder=0, lw=0)


def panel_A(fig, gs, tag="A"):
    set_panel("3A")
    axes = [fig.add_subplot(gs[i]) for i in range(2)]
    for ax, pol in zip(axes, POLS):
        inj = read_table(f"01_qc/{pol}/injections.tsv")
        for _, r in inj.iterrows():
            ax.add_patch(Rectangle((r.injection_order - 0.45, 0), 0.9, 1,
                                   facecolor=ROLE[r.role], edgecolor="#333333", linewidth=0.3))
        spans = block_spans(inj)
        for b, (lo, hi) in spans.items():
            days = pd.to_datetime(inj[inj.block == b].start_time, utc=True).dt.tz_convert(
                "Europe/Vienna")
            lab = f"block {b}: {days.min():%d %b}–{days.max():%d %b %Y}"
            ax.text((lo + hi) / 2, 1.12, lab, ha="center", va="bottom", fontsize=6)
            if lo > 1:
                ax.axvline(lo - 0.5, color="black", lw=0.8)
        ax.set_xlim(0.4, inj.injection_order.max() + 0.6)
        ax.set_ylim(0, 1)
        ax.set_yticks([0.5])
        ax.set_yticklabels([f"{POL_LABEL[pol]} (n = {len(inj)})"])
        ax.tick_params(axis="y", length=0)
        ax.spines["left"].set_visible(False)
    axes[1].set_xlabel("injection order", labelpad=1)
    letter(axes[0], tag, x=-0.07, y=1.55, title="Injection sequence and roles")
    roles = list(ROLE)
    axes[1].legend([Rectangle((0, 0), 1, 1, facecolor=ROLE[k], edgecolor="#333333", lw=0.3)
                    for k in roles], [ROLE_LABEL[k] for k in roles], ncol=6,
                   loc="upper center", bbox_to_anchor=(0.5, -0.9), handlelength=1.0,
                   columnspacing=1.2)


def panel_B(fig, gs, tag="B", compact=False):
    set_panel("3B")
    axes = [fig.add_subplot(gs[i]) for i in range(2)]
    for ax, pol in zip(axes, POLS):
        prov = read_yaml(f"01_qc/{pol}/provenance.yaml")
        thr = prov["params"]["thresholds_from_run"]["qc_feature_rsd_pct"]
        rsd = read_table(f"01_qc/{pol}/rsd_per_feature.tsv")
        summ = read_table(f"01_qc/{pol}/rsd_summary.tsv").set_index("mode")
        sets = [("peak_picked_only", "eligible_peak_picked", "rsd_peak_picked_pct", "-",
                 "peak-picked only"),
                ("including_gap_filled", "eligible_with_gapfill", "rsd_with_gapfill_pct", "--",
                 "incl. gap-filled")]
        for mode, flag, col, ls, lab in sets:
            v = np.sort(rsd.loc[rsd[flag].astype(bool), col].dropna().to_numpy())
            n = int(summ.loc[mode, "n_eligible"])
            if len(v) != n:
                raise ValueError(f"{pol} {mode}: {len(v)} RSD values vs n_eligible {n}")
            y = np.arange(1, len(v) + 1) / len(v)
            sep = "\n" if compact else ": "
            ax.step(v, y, where="post", color=POL[pol], ls=ls, lw=1.0,
                    label=f"{lab}{sep}n = {n:,}, {summ.loc[mode, 'pct_of_eligible_le_30']}% "
                          f"≤ {thr:g}%")
        ax.axvline(thr, color="black", lw=0.6, ls=":")
        ax.set_xlim(0, 150)
        ax.set_ylim(0, 1.0)
        ax.set_xlabel("pooled-QC RSD (%)")
        ax.set_title(POL_LABEL[pol], loc="left", fontsize=7, color=POL[pol], pad=2)
        if compact:
            ax.legend(loc="lower right", handlelength=1.5, fontsize=5, borderaxespad=0.1)
        else:
            ax.legend(loc="center right", bbox_to_anchor=(1.0, 0.22), handlelength=1.8,
                      fontsize=5.8)
    axes[1].tick_params(labelleft=False)
    axes[0].set_ylabel("cumulative fraction" if compact else
                       "cumulative fraction of eligible features")
    letter(axes[0], tag, x=-0.18, y=1.08, title="Pooled-QC RSD per feature")


def panel_C(fig, gs, tag="C"):
    set_panel("3C")
    axes = np.array([[fig.add_subplot(gs[r, c]) for c in range(2)] for r in range(2)])
    for c, pol in enumerate(POLS):
        prov = read_yaml(f"01_qc/{pol}/provenance.yaml")["params"]["thresholds_from_run"]
        rt_tol, ppm_tol = prov["rt_tolerance_min"], prov["mass_accuracy_ppm"]
        inj = read_table(f"01_qc/{pol}/injections.tsv")
        ipi = read_table(f"01_qc/{pol}/is_per_injection.tsv")
        byb = read_table(f"01_qc/{pol}/is_by_block.tsv")
        ipi = ipi[ipi.detected.astype(bool)].merge(
            byb[["compound", "block", "median_rt"]], on=["compound", "block"], how="left")
        ipi["rt_minus_block_median"] = ipi.rt - ipi.median_rt
        wide = byb.pivot(index="compound", columns="block", values="median_rt").dropna()
        shift = (wide[2] - wide[1]).median()
        spans = block_spans(inj)
        n_is = ipi.compound.nunique()
        for r, (col, tol, ylab, lab) in enumerate([
                ("rt_minus_block_median", rt_tol, "RT − block median (min)",
                 f"±{rt_tol:g} min; block 2 − block 1: {shift:+.3f} min (median, n = {len(wide)} IS)"),
                ("ppm_error", ppm_tol, "mass error (ppm)", f"±{ppm_tol:g} ppm")]):
            ax = axes[r, c]
            shade_blocks(ax, spans)
            ax.scatter(ipi.injection_order, ipi[col], s=4, color=POL[pol], alpha=0.7, lw=0,
                       zorder=3)
            for s_ in (-1, 1):
                ax.axhline(s_ * tol, color="black", lw=0.6, ls="--", zorder=2)
            ax.text(0.99, 0.97, lab, transform=ax.transAxes, ha="right", va="top", fontsize=5.8,
                    bbox=dict(facecolor="white", edgecolor="none", pad=0.5, alpha=0.8), zorder=4)
            ax.set_xlim(0.4, inj.injection_order.max() + 0.6)
            lim = max(np.nanmax(np.abs(ipi[col])), tol) * 1.35
            ax.set_ylim(-lim, lim)
            if c == 0:
                ax.set_ylabel(ylab)
            if r == 0:
                ax.set_title(f"{POL_LABEL[pol]} ({n_is} IS detected)", loc="left", fontsize=7,
                             color=POL[pol], pad=2)
                ax.tick_params(labelbottom=False)
            else:
                ax.set_xlabel("injection order")
    letter(axes[0, 0], tag, x=-0.13, y=1.12,
           title="Internal standards: RT and mass accuracy, all injections (shaded = block 1)")


def panel_D(fig, gs, tag="D"):
    set_panel("3D")
    ax = fig.add_subplot(gs)
    data = {pol: read_table(f"01_qc/{pol}/is_vs_calibration.tsv") for pol in POLS}
    if not gates()["reconstitution_volume_confirmed"]:
        pending(ax, "reconstitution\nvolume not confirmed")
        letter(ax, tag, x=0.0, y=1.02, title="IS in soil vs calibration range")
        return "PENDING", "gated on reconstitution_volume_confirmed"
    rows = []
    for pol, d in data.items():
        d = d[d.n_samples_detected > 0].copy()
        d["pol"] = pol
        rows.append(d)
    d = pd.concat(rows).reset_index(drop=True)
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
    letter(ax, tag, x=-0.02, y=1.02, title="IS in soil vs calibration range")
    return "REAL", ""




def rsd_main(fig, gs, tag):
    """Main-text version: peak-picked only (slate) vs including gap-filled values (grey)."""
    from common import SOIL, finding
    set_panel("rsd_main")
    axes = [fig.add_subplot(gs[i]) for i in range(2)]
    stats = {}
    for ax, pol in zip(axes, POLS):
        thr = read_yaml(f"01_qc/{pol}/provenance.yaml")["params"]["thresholds_from_run"][
            "qc_feature_rsd_pct"]
        rsd = read_table(f"01_qc/{pol}/rsd_per_feature.tsv")
        summ = read_table(f"01_qc/{pol}/rsd_summary.tsv").set_index("mode")
        for mode, flag, col, colour, lab in [
                ("peak_picked_only", "eligible_peak_picked", "rsd_peak_picked_pct", SOIL,
                 "detected peaks only"),
                ("including_gap_filled", "eligible_with_gapfill", "rsd_with_gapfill_pct",
                 "#A6A6A6", "incl. gap-filled values")]:
            v = np.sort(rsd.loc[rsd[flag].astype(bool), col].dropna().to_numpy())
            n = int(summ.loc[mode, "n_eligible"])
            if len(v) != n:
                raise ValueError(f"{pol} {mode}: {len(v)} RSD values vs n_eligible {n}")
            pct = summ.loc[mode, "pct_of_eligible_le_30"]
            stats[(pol, mode)] = (n, pct)
            ax.step(v, np.arange(1, n + 1) / n, where="post", color=colour, lw=1.2,
                    label=f"{lab}\n{n:,} features, {pct} % ≤ {thr:g} %")
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
    n_pp, _ = stats[("pos", "peak_picked_only")]
    n_gf, pct_gf = stats[("pos", "including_gap_filled")]
    finding(axes[0], tag, f"Gap filling adds {n_gf / n_pp:.0f}× more features; "
            f"only {pct_gf:.0f} % reproducible (POS)", x=-0.02)
