"""Panels: what the Soil Life Atlas adds to a soil dataset (ClimGrass, positive mode).

Inputs: 06_atlas_match/pos/{coverage,best_match,phylum_detection}.tsv and
02_annotation/pos/{annotations,confidence}.tsv. A soil feature 'matches' when its MS/MS spectrum
matches an Atlas reference spectrum (modified cosine; thresholds in the folder's provenance.yaml);
its occurrence is the set of phyla in which that reference feature was detected.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from common import OI, letter, read_table, read_yaml, set_panel

ATLAS = OI["green"]
ANNOT = OI["blue"]
UNANNOT = "#9AD3BC"  # light green: matched, no annotation


def _matched() -> pd.DataFrame:
    best = read_table("06_atlas_match/pos/best_match.tsv")
    ann = read_table("02_annotation/pos/annotations.tsv", low_memory=False)
    m = best.merge(ann[["feature_id", "confidence", "lipid_class"]],
                   left_on="soil_feature_id", right_on="feature_id", how="left")
    m["annotated"] = m.confidence.notna()
    return m


def coverage(ax) -> str:
    """How far the Atlas reaches into the soil feature table. Returns the panel title."""
    set_panel("atlas_coverage")
    cov = read_table("06_atlas_match/pos/coverage.tsv").iloc[0]
    conf = read_table("02_annotation/pos/confidence.tsv").set_index("confidence")
    m = _matched()
    n_all = int(conf.n_features.sum())
    n_ann = int(conf.loc[["Gold", "Silver", "Bronze"], "n_features"].sum())
    n_ms2 = int(cov.soil_features_with_ms2)
    n_match = int(cov.soil_features_matched)
    if n_match != len(m):
        raise ValueError(f"coverage {n_match} vs best_match rows {len(m)}")
    n_match_ann = int(m.annotated.sum())
    rows = [
        ("features detected", [(n_all, "#CCCCCC", "")]),
        ("with MS/MS spectrum", [(n_ms2, "#AAAAAA", "")]),
        ("matched to the Atlas", [(n_match_ann, ANNOT, "also annotated"),
                                  (n_match - n_match_ann, UNANNOT, "no annotation")]),
        ("annotated (any tool)", [(n_match_ann, ANNOT, "Atlas match"),
                                  (n_ann - n_match_ann, "#CCCCCC", "no match")]),
    ]
    for i, (lab, segs) in enumerate(rows):
        left = 0
        for n, c, _ in segs:
            ax.barh(i, n, left=left, color=c, height=0.62, edgecolor="white", lw=0.4)
            left += n
        total = sum(n for n, _, _ in segs)
        detail = "" if len(segs) == 1 else " (" + ", ".join(
            f"{n:,} {t}" for n, _, t in segs) + ")"
        ax.text(left + n_all * 0.01, i, f"{total:,}{detail}", va="center", fontsize=5.8)
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels([r[0] for r in rows])
    ax.invert_yaxis()
    ax.set_xlim(0, n_all * 1.1)
    ax.set_xlabel("soil features (ClimGrass, positive mode)")
    ax.tick_params(axis="y", length=0)
    return (f"{cov.pct_matched:g}% of soil MS/MS spectra match an Atlas reference; "
            f"{(n_match - n_match_ann) / n_match * 100:.0f}% of the matched features have no "
            "annotation")


def breadth(ax, tag):
    """In how many phyla the matched reference feature occurs."""
    set_panel("atlas_breadth")
    m = _matched()
    n_phyla = int(read_table("06_atlas_match/pos/coverage.tsv").n_phyla_in_presence.iloc[0])
    bins = np.arange(1, n_phyla + 1)
    a = m[m.annotated].n_ref_phyla.value_counts().reindex(bins, fill_value=0)
    u = m[~m.annotated].n_ref_phyla.value_counts().reindex(bins, fill_value=0)
    ax.bar(bins, u, color=UNANNOT, width=0.8, label="no annotation")
    ax.bar(bins, a, bottom=u, color=ANNOT, width=0.8, label="annotated")
    med = m.n_ref_phyla.median()
    ax.axvline(med, color="black", lw=0.6, ls=":")
    ax.text(med + 0.3, ax.get_ylim()[1] * 0.95, f"median {med:g} phyla", fontsize=5.8, va="top")
    ax.set_xticks(bins[::3])
    ax.set_xlim(0.4, n_phyla + 0.6)
    ax.set_xlabel(f"phyla with the matched lipid (of {n_phyla})")
    ax.set_ylabel("matched soil features")
    ax.legend(loc="center right")
    letter(ax, tag, x=-0.25, y=1.02, title="From group-specific to widespread")


def per_phylum(ax, tag):
    """Matched soil features whose reference occurs in each phylum, and only there."""
    set_panel("atlas_per_phylum")
    m = _matched()
    det = read_table("06_atlas_match/pos/phylum_detection.tsv")
    units = det.groupby("phylum").n_units_phylum.first()
    occ = m.ref_phyla.str.split(";").explode().value_counts()
    only = m[m.n_ref_phyla == 1].ref_phyla.value_counts()
    order = occ.sort_values().index
    y = np.arange(len(order))
    ax.barh(y, occ.reindex(order), color=UNANNOT, height=0.7, label="occurs in phylum")
    ax.barh(y, only.reindex(order, fill_value=0), color=ATLAS, height=0.7,
            label="occurs only in this phylum")
    ax.set_yticks(y)
    ax.set_yticklabels([f"{p} ({int(units.get(p, 0))})" for p in order])
    ax.set_xlabel("matched soil features")
    ax.tick_params(axis="y", labelsize=5.8)
    ax.legend(loc="lower right", fontsize=5.5)
    letter(ax, tag, x=-0.55, y=1.02, title="Phyla behind the matches")


def class_breadth(ax, tag, top=10):
    """Annotated matches: lipid class against occurrence breadth."""
    set_panel("atlas_class_breadth")
    m = _matched()
    a = m[m.annotated]
    order = a.lipid_class.value_counts().head(top).index
    order = sorted(order, key=lambda c: a[a.lipid_class == c].n_ref_phyla.median())
    data = [a[a.lipid_class == c].n_ref_phyla.to_numpy() for c in order]
    ax.boxplot(data, vert=False, widths=0.6, showfliers=False,
               medianprops=dict(color="black", lw=0.8), boxprops=dict(lw=0.6),
               whiskerprops=dict(lw=0.6), capprops=dict(lw=0.6))
    rng = np.random.default_rng(0)  # vertical jitter only, for legibility
    for i, d in enumerate(data, 1):
        ax.scatter(d, i + rng.uniform(-0.2, 0.2, len(d)), s=2, color=ANNOT, alpha=0.35, lw=0)
    ax.set_yticks(range(1, len(order) + 1))
    ax.set_yticklabels([f"{c} ({len(d)})" for c, d in zip(order, data)])
    ax.set_xlabel("phyla with the matched lipid")
    letter(ax, tag, x=-0.3, y=1.02, title="Annotated matches by class")


def match_settings() -> str:
    p = read_yaml("06_atlas_match/pos/provenance.yaml")["params"]
    return (f"modified cosine ≥ {p['min_cosine']}, ≥ {p['min_matched_peaks']} shared fragments, "
            f"precursor ± {p['ppm']:g} ppm")
