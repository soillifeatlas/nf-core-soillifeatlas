"""Graphical abstract: the Atlas (reference organisms) -> the soil workflow -> what a soil lipid
profile gains. Every number is read from paper/ (release counts via figures/inputs/)."""
from __future__ import annotations

from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

import panels_atlas as atlas
from common import (ATLAS, GROUP_COLOR, GROUP_ORDER, NAMED, NOMATCH, SOIL, atlas_groups,
                    read_table, set_panel)


def _box(ax, x, y, w, h, fc, ec="none", lw=0.8, ls="-", r=1.5, alpha=1.0):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle=f"round,pad=0,rounding_size={r}",
                                facecolor=fc, edgecolor=ec, lw=lw, linestyle=ls, alpha=alpha))


def _arrow(ax, x0, x1, y, color="#888888"):
    ax.add_patch(FancyArrowPatch((x0, y), (x1, y), arrowstyle="simple,head_width=4,head_length=4",
                                 mutation_scale=2.2, color=color, lw=0))


def graphical_abstract(ax):
    set_panel("ga")
    ax.set_xlim(0, 200)
    ax.set_ylim(0, 80)
    ax.set_axis_off()

    # ------------------------------------------------------------------ headline
    ax.text(2, 77, "Which organisms could have made the lipids in a soil?", fontsize=10,
            fontweight="bold", va="top", color=SOIL)
    ax.text(2, 70.5, "The Soil Life Atlas links soil lipid spectra to reference organisms; an "
            "open workflow measures, checks and annotates the soil data", fontsize=6.8,
            va="top", color="#444444")

    # ------------------------------------------------------------------ left: the Atlas
    g = atlas_groups()
    g = g[g.analysed.astype(str) == "True"]
    per_group = g.groupby("group").samples_pos.sum().reindex(GROUP_ORDER)
    n_phyla = len(g)
    _box(ax, 2, 4, 58, 58, ATLAS, alpha=0.10, ec=ATLAS, lw=0.9)
    ax.text(6, 58.5, "Soil Life Atlas", fontsize=8.5, fontweight="bold", color=ATLAS, va="top")
    ax.text(6, 53.5, f"{int(per_group.sum())} reference samples in {n_phyla} phyla\n"
            "spectra · response · occurrence", fontsize=5.6, color="#335", va="top",
            linespacing=1.3)
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

    # ------------------------------------------------------------------ middle: the workflow
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
        ax.text(75, y + 1.4, t, fontsize=6.4, fontweight="bold", va="center",
                color=SOIL if nf else ATLAS)
        ax.text(75, y - 1.9, d.replace("\n", " "), fontsize=5.2, va="center", color="#555")
        if i:
            ax.annotate("", xy=(94, y + 4.1), xytext=(94, y + 6.1),
                        arrowprops=dict(arrowstyle="-|>", lw=0.6, color="#888",
                                        mutation_scale=5))
    ax.text(72, 7.5, "open nf-core workflow · laptop or HPC", fontsize=5.4, color="#666")

    _arrow(ax, 118, 126.5, 33, SOIL)

    # ------------------------------------------------------------------ right: what is gained
    set_panel("ga")
    cov = read_table("06_atlas_match/pos/coverage.tsv").iloc[0]
    m = atlas._matched()
    n_ms2, n_match = int(cov.soil_features_with_ms2), int(cov.soil_features_matched)
    n_named = int(m.named.sum())
    ax.text(128, 58.5, "What the soil profile gains", fontsize=8.5, fontweight="bold",
            color=SOIL, va="top")
    ax.text(128, 53.5, "ClimGrass grassland soils, positive mode", fontsize=5.6, color="#444",
            va="top")
    x, W = 128, 68
    ax.text(x, 49.5, f"{cov.pct_matched:.0f} %", fontsize=14, fontweight="bold", color=ATLAS,
            va="top")
    ax.text(x + 18, 47, f"of soil MS/MS spectra match an Atlas\nreference ({n_match:,} of "
            f"{n_ms2:,})", fontsize=5.6, va="center", linespacing=1.2)
    segs = [(n_match - n_named, ATLAS), (n_named, NAMED), (n_ms2 - n_match, NOMATCH)]
    left = x
    for n, c in segs:
        w = W * n / n_ms2
        ax.add_patch(FancyBboxPatch((left, 37.5), w, 3.6, boxstyle="square,pad=0", facecolor=c,
                                    edgecolor="white", lw=0.6))
        left += w
    ax.text(x, 35.5, f"{(n_match - n_named) / n_match * 100:.0f} % of the matched lipids have "
            "no name;\nthe Atlas still shows where they occur", fontsize=5.6, va="top",
            color="#0B5F57", linespacing=1.2)
    occ = m.explode("groups").groups.value_counts().reindex(GROUP_ORDER, fill_value=0)
    ax.text(x, 27, "matched soil lipids found in each group", fontsize=5.8, va="top",
            color="#444")
    for i, grp in enumerate(GROUP_ORDER):
        y = 21.5 - i * 3.9
        pct = occ[grp] / len(m) * 100
        ax.text(x, y, grp, fontsize=5.6, va="center")
        _box(ax, x + 17, y - 1.2, 42 * pct / 100, 2.4, GROUP_COLOR[grp], r=0.6, alpha=0.9)
        ax.text(x + 18 + 42 * pct / 100, y, f"{pct:.0f} %", fontsize=5.4, va="center")
