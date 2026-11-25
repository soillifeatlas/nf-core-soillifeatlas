"""Panels: interpreting individual soil features with the Atlas (mirror plots, occurrence)."""
from __future__ import annotations

import numpy as np

from common import OI, POL, letter, read_table, set_panel

ATLAS = OI["green"]
EX_COLOR = {"broad": OI["blue"], "restricted": OI["vermillion"]}


def _phyla(n: int) -> str:
    return f"{n} phylum" if n == 1 else f"{n} phyla"


def mirror(ax, role: str, tag: str):
    """Soil MS/MS spectrum (top) against its Atlas reference spectrum (bottom)."""
    set_panel(f"mirror_{role}")
    exs = read_table("06_atlas_match/examples/examples.tsv")
    ex = exs[exs.role == role].iloc[0]
    m = read_table(f"06_atlas_match/examples/mirror_{ex.soil_feature_id}.tsv")
    soil, atlas = m[m.side == "soil"], m[m.side == "atlas"]
    # rel_intensity is already in % of the base peak
    ax.vlines(soil.mz, 0, soil.rel_intensity, color=POL["pos"], lw=0.6)
    ax.vlines(atlas.mz, 0, -atlas.rel_intensity, color=ATLAS, lw=0.6)
    ax.axhline(0, color="black", lw=0.5)
    ax.set_ylim(-128, 128)
    ax.set_yticks([-100, -50, 0, 50, 100])
    ax.set_yticklabels(["100", "50", "0", "50", "100"])
    ax.set_xlabel("m/z")
    ax.set_ylabel("relative intensity (%)")
    ax.text(0.01, 0.98, f"soil: m/z {ex.soil_mz:.4f}, RT {ex.rt:.2f} min", transform=ax.transAxes,
            va="top", fontsize=6, color=POL["pos"])
    ax.text(0.01, 0.02, f"Atlas reference: m/z {ex.ref_mz:.4f}, found in "
            f"{_phyla(int(ex.n_ref_phyla))}", transform=ax.transAxes, va="bottom", fontsize=6,
            color=ATLAS)
    ax.text(0.99, 0.56, f"cosine {ex.cosine:.3f}\n{int(ex.matched_peaks)} shared fragments",
            transform=ax.transAxes, va="bottom", ha="right", fontsize=6,
            bbox=dict(facecolor="white", edgecolor="none", pad=0.6, alpha=0.85))
    letter(ax, tag, x=-0.1, y=1.03,
           title=f"{ex.lipid_name}: soil spectrum vs Atlas reference ({role} occurrence)")
    return ex


def phylum_detection(ax, tag: str):
    set_panel("phylum_detection")
    exs = read_table("06_atlas_match/examples/examples.tsv")
    det = read_table("06_atlas_match/examples/examples_phylum_detection.tsv")
    units = det.groupby("phylum").n_units_phylum.first()
    broad = exs[exs.role == "broad"].ref_feature_id.iloc[0]
    order = (det[det.feature_id == broad].sort_values("detection_rate").phylum.tolist())
    order = [p for p in units.index if p not in order] + order
    y = np.arange(len(order))
    for k, (_, ex) in enumerate(exs.iterrows()):
        off = 0.2 if k == 0 else -0.2
        d = det[det.feature_id == ex.ref_feature_id].set_index("phylum").reindex(order)
        rate = d.detection_rate.fillna(0).to_numpy() * 100
        ax.barh(y + off, rate, height=0.4, color=EX_COLOR[ex.role],
                label=f"{ex.lipid_name} ({_phyla(int(ex.n_ref_phyla))})")
        for yi, (rt, nd) in enumerate(zip(rate, d.n_detected)):
            if nd == nd and nd > 0:
                ax.text(rt + 1, yi + off, f"{int(nd)}/{int(units[order[yi]])}", va="center",
                        fontsize=4.8)
    ax.set_yticks(y)
    ax.set_yticklabels(order)
    ax.set_xlim(0, 118)
    ax.set_xlabel("reference samples with the lipid (%)")
    ax.legend(loc="lower right", fontsize=5.8)
    letter(ax, tag, x=-0.32, y=1.03, title="Where the two lipids occur in the Atlas")
