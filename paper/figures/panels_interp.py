"""Panels: interpreting individual soil features with the Atlas (mirror plots, occurrence)."""
from __future__ import annotations

import numpy as np

from common import (ATLAS, GROUP_COLOR, GROUP_ORDER, SOIL, atlas_groups, finding, read_table,
                    set_panel)

FRAG_TOL_DA = 0.02  # same fragment tolerance as the matching (06_atlas_match provenance)


def _phyla(n: int) -> str:
    return f"{n} phylum" if n == 1 else f"{n} phyla"


def _example(role: str):
    exs = read_table("06_atlas_match/examples/examples.tsv")
    return exs[exs.role == role].iloc[0]


def _shared_fragments(soil, atlas, n=3):
    """Top-n fragments present in both spectra (within FRAG_TOL_DA), ranked by the weaker side."""
    out = []
    a_mz = atlas.mz.to_numpy()
    for _, r in soil.iterrows():
        j = np.argmin(np.abs(a_mz - r.mz))
        if abs(a_mz[j] - r.mz) <= FRAG_TOL_DA:
            out.append((min(r.rel_intensity, atlas.rel_intensity.iloc[j]), r.mz,
                        r.rel_intensity))
    return sorted(out, reverse=True)[:n]


def mirror(ax, role: str, tag: str, headline: str):
    """Soil MS/MS spectrum (up) against its Atlas reference spectrum (down), zoomed to the
    precursor range; the strongest shared fragments are labelled with their m/z."""
    set_panel(f"mirror_{role}")
    ex = _example(role)
    m = read_table(f"06_atlas_match/examples/mirror_{ex.soil_feature_id}.tsv")
    soil, atlas = m[m.side == "soil"], m[m.side == "atlas"]
    xmax = ex.soil_mz + 25
    hidden = int((m.mz > xmax).sum())
    # rel_intensity is already in % of the base peak
    ax.vlines(soil.mz, 0, soil.rel_intensity, color=SOIL, lw=0.7)
    ax.vlines(atlas.mz, 0, -atlas.rel_intensity, color=ATLAS, lw=0.7)
    ax.axhline(0, color="black", lw=0.5)
    placed = []
    for _, mz, h in _shared_fragments(soil, atlas, n=4):
        if any(abs(mz - p) < 45 for p in placed):  # keep labels from overlapping
            continue
        placed.append(mz)
        ax.text(mz, h + 4, f"{mz:.2f}", ha="center", va="bottom", fontsize=5.5)
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
    if hidden:
        note += f"\n{hidden} peak{'s' if hidden > 1 else ''} above m/z {xmax:.0f} not shown"
    ax.text(0.99, 0.36, note, transform=ax.transAxes, va="top", ha="right", fontsize=5.8,
            color="#333333")
    finding(ax, tag, headline.format(name=ex.lipid_name, phyla=_phyla(int(ex.n_ref_phyla))),
            x=-0.02)
    return ex


def occurrence_matrix(ax, tag: str):
    """Both example lipids across every analysed phylum, phyla grouped by organism group."""
    set_panel("phylum_detection")
    exs = read_table("06_atlas_match/examples/examples.tsv")
    det = read_table("06_atlas_match/examples/examples_phylum_detection.tsv")
    units_all = read_table("06_atlas_match/pos/phylum_detection.tsv") \
        .groupby("phylum").n_units_phylum.first()
    grp = atlas_groups().set_index("phylum").group
    cols = [(g, p) for g in GROUP_ORDER for p in grp[grp == g].index if p in units_all.index]
    exs = exs.set_index("role").loc[["broad", "restricted"]].reset_index()
    for i, (_, ex) in enumerate(exs.iterrows()):
        d = det[det.feature_id == ex.ref_feature_id].set_index("phylum")
        for j, (g, p) in enumerate(cols):
            if p in d.index:
                r = d.loc[p]
                ax.scatter(j, i, s=20 + r.detection_rate * 130, color=GROUP_COLOR[g], lw=0,
                           zorder=3)
                ax.text(j, i + 0.36, f"{int(r.n_detected)}/{int(r.n_units_phylum)}",
                        ha="center", va="top", fontsize=4.6, color="#444444")
            else:
                ax.scatter(j, i, s=9, facecolor="white", edgecolor="#BBBBBB", lw=0.5, zorder=3)
    # group bands
    start = 0
    for g in GROUP_ORDER:
        n = sum(1 for gg, _ in cols if gg == g)
        if not n:
            continue
        ax.add_patch(plt_rect(start - 0.45, -0.95, n - 0.1, 0.22, GROUP_COLOR[g]))
        ax.text(start + (n - 1) / 2, -1.02, g, ha="center", va="bottom", fontsize=6,
                fontweight="bold", color=GROUP_COLOR[g])
        start += n
    ax.set_xticks(range(len(cols)))
    ax.set_xticklabels([p for _, p in cols], rotation=45, ha="right", fontsize=5.8)
    ax.set_yticks(range(len(exs)))
    ax.set_yticklabels([f"{r.lipid_name}\n({_phyla(int(r.n_ref_phyla))})"
                        for _, r in exs.iterrows()])
    ax.set_xlim(-0.6, len(cols) - 0.4)
    ax.set_ylim(len(exs) - 0.35, -1.25)
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
    for _, ex in exs.iterrows():
        gs = sorted({grp[p] for p in det[det.feature_id == ex.ref_feature_id].phylum})
        where = gs[0].lower() if len(gs) == 1 else f"{len(gs)} of {len(GROUP_ORDER)} groups"
        parts.append(f"{ex.lipid_name} occurs in {where}")
    finding(ax, tag, "; ".join(parts), x=-0.02, y=1.02)


def plt_rect(x, y, w, h, color):
    from matplotlib.patches import Rectangle
    return Rectangle((x, y), w, h, facecolor=color, alpha=0.25, edgecolor="none")
