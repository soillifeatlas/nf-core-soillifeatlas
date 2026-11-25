"""Panels: what the Soil Life Atlas adds to a soil dataset (ClimGrass, positive mode).

Inputs: 06_atlas_match/pos/{coverage,best_match,phylum_detection,provenance}, 02_annotation/pos/
annotations.tsv and the phylum -> organism-group mapping (figures/inputs/atlas_groups.csv).
A soil feature 'matches' when its MS/MS spectrum matches an Atlas reference spectrum; it
'occurs' in a group when the matched reference feature was detected in a phylum of that group.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from common import (ATLAS, GROUP_COLOR, GROUP_ORDER, NAMED, NOMATCH, atlas_groups, finding,
                    read_table, read_yaml, set_panel)


def _matched() -> pd.DataFrame:
    best = read_table("06_atlas_match/pos/best_match.tsv")
    ann = read_table("02_annotation/pos/annotations.tsv", low_memory=False)
    grp = atlas_groups().set_index("phylum").group
    m = best.merge(ann[["feature_id", "confidence", "lipid_class"]],
                   left_on="soil_feature_id", right_on="feature_id", how="left")
    m["named"] = m.confidence.notna()
    phyla = m.ref_phyla.str.split(";")
    missing = set(phyla.explode()) - set(grp.index)
    if missing:
        raise ValueError(f"phyla without an organism group: {missing}")
    m["groups"] = phyla.apply(lambda ps: sorted({grp[p] for p in ps}))
    m["n_groups"] = m.groups.str.len()
    return m


def _group_samples() -> pd.Series:
    """Reference samples per group in the matching presence table (matching denominators)."""
    det = read_table("06_atlas_match/pos/phylum_detection.tsv")
    grp = atlas_groups().set_index("phylum").group
    return det.groupby("phylum").n_units_phylum.first().groupby(grp).sum()


def coverage(ax, tag):
    set_panel("atlas_coverage")
    cov = read_table("06_atlas_match/pos/coverage.tsv").iloc[0]
    conf = read_table("02_annotation/pos/confidence.tsv")
    m = _matched()
    n_all = int(conf.n_features.sum())
    n_ms2, n_match = int(cov.soil_features_with_ms2), int(cov.soil_features_matched)
    if n_match != len(m):
        raise ValueError(f"coverage {n_match} vs best_match rows {len(m)}")
    n_named = int(m.named.sum())
    segs = [(n_match - n_named, ATLAS, "matched, no name", "white"),
            (n_named, NAMED, "matched and named", "white"),
            (n_ms2 - n_match, NOMATCH, "no Atlas match", "#333333")]
    left = 0
    for n, c, lab, tc in segs:
        ax.barh(0, n, left=left, color=c, height=0.6, edgecolor="white", lw=0.8)
        ax.text(left + n / 2, 0, f"{n:,}", ha="center", va="center", color=tc, fontsize=7,
                fontweight="bold")
        ax.text(left + n / 2, 0.45, f"{lab}\n{n / n_ms2 * 100:.0f} %", ha="center",
                va="bottom", fontsize=6)
        left += n
    ax.set_xlim(0, n_ms2)
    ax.set_ylim(-0.45, 1.25)
    ax.set_axis_off()
    ax.text(0, -0.42, f"{n_ms2:,} soil features with an MS/MS spectrum (of {n_all:,} detected); "
            f"ClimGrass, positive mode", fontsize=6, color="#555555", va="top")
    finding(ax, tag, f"{cov.pct_matched:.0f} % of soil spectra match an Atlas reference; "
            f"{(n_match - n_named) / n_match * 100:.0f} % of these have no name",
            x=0.0, y=1.08)


def breadth(ax, tag):
    set_panel("atlas_breadth")
    m = _matched()
    bins = np.arange(1, len(GROUP_ORDER) + 1)
    u = m[~m.named].n_groups.value_counts().reindex(bins, fill_value=0)
    n = m[m.named].n_groups.value_counts().reindex(bins, fill_value=0)
    ax.bar(bins, u, color=ATLAS, width=0.75, label="no name")
    ax.bar(bins, n, bottom=u, color=NAMED, width=0.75, label="named")
    one = (m.n_groups == 1).mean() * 100
    ax.set_xticks(bins)
    ax.set_xlabel("organism groups in which the lipid occurs")
    ax.set_ylabel("matched soil features")
    ax.legend(loc="upper right", handlelength=1)
    finding(ax, tag, f"{one:.0f} % point to a single organism group", x=-0.02)


def per_group(ax, tag):
    set_panel("atlas_per_group")
    m = _matched()
    samples = _group_samples()
    occ = m.explode("groups").groups.value_counts().reindex(GROUP_ORDER, fill_value=0)
    only = m[m.n_groups == 1].groups.str[0].value_counts().reindex(GROUP_ORDER, fill_value=0)
    order = occ.sort_values().index
    y = np.arange(len(order))
    for yi, g in zip(y, order):
        ax.barh(yi, occ[g] / len(m) * 100, color=GROUP_COLOR[g], alpha=0.35, height=0.7)
        ax.barh(yi, only[g] / len(m) * 100, color=GROUP_COLOR[g], height=0.7)
    ax.set_yticks(y)
    ax.set_yticklabels([f"{g} ({int(samples[g])})" for g in order])
    ax.tick_params(axis="y", length=0)
    ax.set_xlabel("matched soil features (%)")
    ax.text(0.98, 0.04, "pale: occurs in group\nsolid: only in this group\n(n reference samples)",
            transform=ax.transAxes, ha="right", va="bottom", fontsize=5.5, color="#444444")
    finding(ax, tag, "Matches point to all six organism groups", x=-0.02)


def class_by_group(ax, tag, top=8):
    set_panel("atlas_class_group")
    m = _matched()
    a = m[m.named]
    order = a.lipid_class.value_counts().head(top).index.tolist()
    e = a[a.lipid_class.isin(order)].explode("groups")
    share = (e.groupby(["lipid_class", "groups"]).soil_feature_id.nunique().unstack(fill_value=0)
             .reindex(index=order, columns=GROUP_ORDER, fill_value=0)
             .div(a.lipid_class.value_counts()[order], axis=0) * 100)
    for i, c in enumerate(order):
        for j, g in enumerate(GROUP_ORDER):
            ax.scatter(j, i, s=share.loc[c, g] * 1.1, color=GROUP_COLOR[g], lw=0)
    ax.set_xticks(range(len(GROUP_ORDER)))
    ax.set_xticklabels(GROUP_ORDER, rotation=35, ha="right")
    ax.set_yticks(range(len(order)))
    ax.set_yticklabels([f"{c} ({int((a.lipid_class == c).sum())})" for c in order])
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
    p = read_yaml("06_atlas_match/pos/provenance.yaml")["params"]
    return (f"modified cosine ≥ {p['min_cosine']}, ≥ {p['min_matched_peaks']} shared fragments, "
            f"precursor ± {p['ppm']:g} ppm")
