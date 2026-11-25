"""Figure 4 - annotation (ClimGrass)."""
from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np

from common import (CONF, DOUBLE, MM, POL, POL_LABEL, Panels, letter, read_table, save,
                    set_panel)

POLS = ("pos", "neg")
SOURCE_ORDER = ["LipidSearch", "SIRIUS/CANOPUS", "MZmine lipid rules", "spectral library"]


def upset(fig, gs, pol, first, tag="A"):
    sets = read_table(f"02_annotation/{pol}/upset_sets.tsv")
    inter = read_table(f"02_annotation/{pol}/upset_intersections.tsv")
    inter = inter.sort_values("n_features", ascending=False).reset_index(drop=True)
    names = [s for s in SOURCE_ORDER if s in set(sets.source)] + \
        [s for s in sets.source if s not in SOURCE_ORDER]
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
    ax_bar.set_title(f"{POL_LABEL[pol]}: {inter.n_features.sum():,} features with a source",
                     loc="left", fontsize=7, color=POL[pol], pad=2)
    ypos = {n: i for i, n in enumerate(names)}
    for xi, combo in zip(x, inter.intersection):
        members = [m.strip() for m in combo.split(";")]
        ys = [ypos[m] for m in members]
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
               title="Annotation sources (UpSet; exclusive intersections; left: set sizes)")


def panel_B(ax, tag="B"):
    set_panel("4B")
    ys = []
    for i, pol in enumerate(POLS):
        c = read_table(f"02_annotation/{pol}/confidence.tsv").set_index("confidence")
        left = 0
        for lvl in ("Gold", "Silver", "Bronze"):
            n = int(c.loc[lvl, "n_features"])
            ax.barh(i, n, left=left, color=CONF[lvl], height=0.6, edgecolor="white", lw=0.4)
            ax.text(left + n / 2, i, f"{n:,}", ha="center", va="center", fontsize=5.5,
                    color="black" if lvl != "Bronze" else "white")
            left += n
        total_all = int(c.n_features.sum())
        ax.text(left, i, f"  {left:,} of {total_all:,}", va="center", fontsize=5.5)
        ys.append(i)
    ax.set_yticks(ys)
    ax.set_yticklabels([POL_LABEL[p] for p in POLS])
    for t, p in zip(ax.get_yticklabels(), POLS):
        t.set_color(POL[p])
    ax.invert_yaxis()
    ax.set_xlabel("annotated features")
    xmax = ax.get_xlim()[1]
    ax.set_xlim(0, xmax * 1.3)
    handles = [plt.Rectangle((0, 0), 1, 1, color=CONF[k]) for k in ("Gold", "Silver", "Bronze")]
    ax.legend(handles, ["Gold", "Silver", "Bronze"], ncol=3, loc="lower right",
              bbox_to_anchor=(1.0, 1.0), handlelength=1, columnspacing=0.8)
    letter(ax, tag, x=-0.15, y=1.15, title="Confidence (ClimGrass)")


def panel_C(fig, gs, top=12, tag="C"):
    set_panel("4C")
    axes = [fig.add_subplot(gs[i]) for i in range(2)]
    for ax, pol in zip(axes, POLS):
        c = read_table(f"02_annotation/{pol}/classes.tsv")
        c = c.sort_values("n_features", ascending=False).head(top)
        x = np.arange(len(c))
        ax.bar(x - 0.2, c.n_features, width=0.4, color=POL[pol], label="features")
        ax.bar(x + 0.2, c.n_unique_lipids, width=0.4, color="white", edgecolor=POL[pol],
               lw=0.7, label="unique lipids")
        ax.set_xticks(x)
        ax.set_xticklabels(c.lipid_class, rotation=90)
        ax.set_title(f"{POL_LABEL[pol]}: {top} largest classes", loc="left", fontsize=7,
                     color=POL[pol], pad=2)
        ax.legend(loc="upper right", handlelength=1)
    axes[0].set_ylabel("Gold + Silver + Bronze")
    letter(axes[0], tag, x=-0.12, y=1.1, title="Classes: features vs unique lipids")


def panel_D(fig, gs, tag="D"):
    set_panel("4D")
    ax = fig.add_subplot(gs)
    lo, hi = np.inf, -np.inf
    notes = []
    for pol in POLS:
        d = read_table(f"02_annotation/{pol}/rt_model_PC.tsv")
        t = d[d.rt_check.isin(["ok", "outlier"])]
        ok, out = t[t.rt_check == "ok"], t[t.rt_check == "outlier"]
        ax.scatter(ok.rt_expected, ok.rt, s=5, color=POL[pol], lw=0, alpha=0.8,
                   label=f"{POL_LABEL[pol]} within tolerance ({len(ok)})")
        ax.scatter(out.rt_expected, out.rt, s=9, marker="x", color=POL[pol], lw=0.6,
                   label=f"{POL_LABEL[pol]} outlier ({len(out)})")
        lo, hi = min(lo, t.rt.min(), t.rt_expected.min()), max(hi, t.rt.max(), t.rt_expected.max())
        notes.append(f"{POL_LABEL[pol]} {int((d.rt_check == 'untested').sum())}")
    ax.plot([lo, hi], [lo, hi], color="black", lw=0.6, ls="--")
    ax.set_xlabel("fitted RT (min)")
    ax.set_ylabel("observed RT (min)")
    ax.legend(loc="upper left", fontsize=5.5, handletextpad=0.2, borderaxespad=0.1)
    ax.text(0.02, 0.5, "untested, not drawn:\n" + ", ".join(notes), transform=ax.transAxes, ha="left",
            fontsize=5.5)
    letter(ax, tag, x=-0.15, y=1.04, title="PC retention-time model")


