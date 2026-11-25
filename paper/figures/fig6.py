"""Figure 6 - soil-to-Atlas matching (POS)."""
from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np

from common import (DOUBLE, MM, OI, POL, Panels, letter, pending, read_table, save, set_panel)

ATLAS = OI["green"]
EX_COLOR = {"broad": OI["blue"], "restricted": OI["vermillion"]}


def mirror(ax, ex, tag, title):
    m = read_table(f"06_atlas_match/examples/mirror_{ex.soil_feature_id}.tsv")
    soil, atlas = m[m.side == "soil"], m[m.side == "atlas"]
    ax.vlines(soil.mz, 0, soil.rel_intensity, color=POL["pos"], lw=0.6)
    ax.vlines(atlas.mz, 0, -atlas.rel_intensity, color=ATLAS, lw=0.6)
    ax.axhline(0, color="black", lw=0.5)
    ax.set_ylim(-128, 128)  # rel_intensity is already in % of base peak
    ax.set_yticks([-100, -50, 0, 50, 100])
    ax.set_yticklabels(["100", "50", "0", "50", "100"])
    ax.set_xlabel("m/z")
    ax.set_ylabel("relative intensity (%)")
    ax.text(0.01, 0.98, f"soil feature {ex.soil_feature_id} (m/z {ex.soil_mz:.4f}, "
            f"RT {ex.rt:.2f} min, {ex.confidence})", transform=ax.transAxes, va="top",
            fontsize=5.8, color=POL["pos"])
    ax.text(0.01, 0.02, f"Atlas reference feature {ex.ref_feature_id} (m/z {ex.ref_mz:.4f}, "
            f"{ex.ppm:+.1f} ppm)", transform=ax.transAxes, va="bottom", fontsize=5.8,
            color=ATLAS)
    ax.text(0.99, 0.56, f"cosine {ex.cosine:.3f}\n{int(ex.matched_peaks)} matched fragments",
            transform=ax.transAxes, va="bottom", ha="right", fontsize=6,
            bbox=dict(facecolor="white", edgecolor="none", pad=0.6, alpha=0.85))
    letter(ax, tag, x=-0.1, y=1.03, title=title)


def build():
    P = Panels(6)
    fig = plt.figure(figsize=(DOUBLE, 150 * MM))
    gs = fig.add_gridspec(2, 2, width_ratios=[1.25, 1], height_ratios=[1, 1], hspace=0.45,
                          wspace=0.3, left=0.08, right=0.98, top=0.94, bottom=0.08)
    for tag, role, r in (("A", "broad", 0), ("B", "restricted", 1)):
        set_panel(f"6{tag}")
        exs = read_table("06_atlas_match/examples/examples.tsv")
        ex = exs[exs.role == role].iloc[0]
        mirror(fig.add_subplot(gs[r, 0]), ex, tag,
               f"{ex.lipid_name}: soil (top) vs Atlas reference (bottom), {role} example")
        P.add(tag, "REAL", f"{ex.lipid_name}; examples proposed, to confirm")

    set_panel("6C")
    ax = fig.add_subplot(gs[0, 1])
    exs = read_table("06_atlas_match/examples/examples.tsv")
    det = read_table("06_atlas_match/examples/examples_phylum_detection.tsv")
    units = det.groupby("phylum").n_units_phylum.first()
    phyla = units.index.tolist()
    order = (det[det.feature_id == exs[exs.role == "broad"].ref_feature_id.iloc[0]]
             .sort_values("detection_rate", ascending=True).phylum.tolist())
    order = [p for p in phyla if p not in order] + order
    y = np.arange(len(order))
    for k, (_, ex) in enumerate(exs.iterrows()):
        d = det[det.feature_id == ex.ref_feature_id].set_index("phylum").reindex(order)
        rate = d.detection_rate.fillna(0).to_numpy() * 100
        ax.barh(y + (0.2 if k == 0 else -0.2), rate, height=0.4, color=EX_COLOR[ex.role],
                label=f"{ex.lipid_name} ({int(ex.n_ref_phyla)} "
                      f"{'phylum' if int(ex.n_ref_phyla) == 1 else 'phyla'})")
        for yi, (rt, nd) in enumerate(zip(rate, d.n_detected)):
            if nd == nd and nd > 0:
                ax.text(rt + 1, yi + (0.2 if k == 0 else -0.2),
                        f"{int(nd)}/{int(units[order[yi]])}", va="center", fontsize=4.8)
    ax.set_yticks(y)
    ax.set_yticklabels(order)
    ax.set_xlim(0, 118)
    ax.set_xlabel("reference samples with detection (%)")
    ax.legend(loc="lower right", fontsize=5.8)
    letter(ax, "C", x=-0.32, y=1.03, title="Detection per phylum (Atlas release reprocessing)")
    P.add("C", "REAL", "only phyla where either reference feature is detected are listed")

    set_panel("6D")
    axD = fig.add_subplot(gs[1, 1])
    pending(axD, "Explorer screenshot\n(issue #7)")
    letter(axD, "D", x=0.0, y=1.03, title="Atlas Explorer record")
    P.add("D", "PENDING", "Explorer screenshot, issue #7", sources=[])
    return save(fig, 6, "interpretation"), P.rows
