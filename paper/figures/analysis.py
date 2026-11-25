#!/usr/bin/env python
"""Every number and table behind the manuscript figures and captions.

    python paper/figures/analysis.py            # write figures/derived/ and numbers.tsv
    python paper/figures/analysis.py --check    # also compare with expected_numbers.csv

This is the only place where figure numbers are made. It reads the exported pipeline tables
under paper/ (and the Atlas group mapping in figures/inputs/) and writes:

  figures/derived/<table>.tsv   plot-ready tables, one per panel input
  figures/derived/numbers.tsv   every derived number used in a panel title or caption
  figures/derived/sources.yaml  which paper/ files each table and number came from

plots.py and captions.py only read figures/derived/. Sections follow the figures.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

sys.path.insert(0, str(Path(__file__).resolve().parent))

from common import (FIG_DIR, GROUP_ORDER, atlas_groups, gates, read_table, read_yaml,  # noqa: E402
                    reads, reset_reads, set_panel)

DERIVED = FIG_DIR / "derived"
EXPECTED = FIG_DIR / "expected_numbers.csv"
POLS = ("pos", "neg")
FRAG_TOL_DA = 0.02  # fragment tolerance of the spectral matching (06_atlas_match provenance)

TABLES: dict[str, pd.DataFrame] = {}
NUMBERS: dict[str, object] = {}
SOURCES: dict[str, list[str]] = {}


def table(name: str, df: pd.DataFrame) -> None:
    """Record a plot-ready table with the paper/ files read so far in the current section."""
    import common
    TABLES[name] = df.reset_index(drop=True)
    SOURCES[name] = sorted(set(reads().get(name, [])) |
                           set(reads().get(common._current_panel[0], [])))


def number(key: str, value, group: str) -> None:
    """Record a number; `group` is the set_panel() bucket whose reads it came from."""
    if isinstance(value, (np.integer,)):
        value = int(value)
    if isinstance(value, (np.floating,)):
        value = float(value)
    NUMBERS[key] = value
    SOURCES[f"number:{key}"] = reads().get(group, [])


# =========================================================================== the Atlas (Fig 1)
LAYER_CHAIN = ["consensus", "quality_filtered", "enriched", "with_msms", "annotated"]


def release():
    """Fig 1: what the Atlas release holds (exported by export_groups.py into inputs/)."""
    set_panel("release")
    prov = read_yaml("figures/inputs/atlas_groups_provenance.yaml")
    number("release_label", prov["release"]["label"], "release")
    g = atlas_groups()
    g["analysed"] = g.analysed.astype(str) == "True"
    g["group_rank"] = g.group.map({x: i for i, x in enumerate(GROUP_ORDER)})
    table("release_phyla", g.sort_values(["group_rank", "samples_pos"], ascending=[True, False])
          [["group", "phylum", "samples_pos", "samples_neg", "analysed"]])
    c = read_table("figures/inputs/atlas_release_counts.csv").set_index("key").value
    for k in ("reference_samples_pos", "reference_samples_neg", "phyla_collected",
              "phyla_analysed", "genera", "groups", "lipids_pos", "lipids_neg"):
        number(f"release_{k}", int(c[k]), "release")
    assert int(c.reference_samples_pos) == int(g.samples_pos.sum()), "release counts disagree"
    assert int(c.phyla_collected) == len(g) and int(c.phyla_analysed) == int(g.analysed.sum())
    excl, considered = int(c.key_finding_exclusive_features), int(c.key_finding_features_considered)
    number("release_exclusive_features", excl, "release")
    number("release_features_considered", considered, "release")
    number("release_exclusive_pct", round(excl / considered * 100, 1), "release")
    lay = read_table("figures/inputs/atlas_feature_layers.csv").set_index("id").loc[LAYER_CHAIN]
    table("release_layers", lay.reset_index()[["id", "label", "pos", "neg"]])


# =========================================================================== quality (S1, Fig 2B)
def qc():
    set_panel("qc_sequence")
    seq, blocks = [], []
    for pol in POLS:
        inj = read_table(f"01_qc/{pol}/injections.tsv")
        seq.append(inj[["injection_order", "role", "block"]].assign(pol=pol))
        for b, sub in inj.groupby("block"):
            days = pd.to_datetime(sub.start_time, utc=True).dt.tz_convert("Europe/Vienna")
            blocks.append({"pol": pol, "block": b, "lo": int(sub.injection_order.min()),
                           "hi": int(sub.injection_order.max()),
                           "label": f"block {b}: {days.min():%d %b}–{days.max():%d %b %Y}"})
    table("qc_sequence", pd.concat(seq))
    table("qc_blocks", pd.DataFrame(blocks))
    inj0 = read_table("01_qc/pos/injections.tsv")
    t = pd.to_datetime(inj0.start_time, utc=True)
    number("qc_block_gap_h_pos",
           (t[inj0.block == 2].min() - t[inj0.block == 1].max()).total_seconds() / 3600,
           "qc_sequence")
    number("qc_n_pooled_qc_pos", int((inj0.role == "pooled_qc").sum()), "qc_sequence")
    for pol in POLS:
        number(f"qc_n_injections_{pol}", len(read_table(f"01_qc/{pol}/injections.tsv")),
               "qc_sequence")

    set_panel("qc_rsd")
    vals, summ = [], []
    for pol in POLS:
        thr = read_yaml(f"01_qc/{pol}/provenance.yaml")["params"]["thresholds_from_run"]
        rsd = read_table(f"01_qc/{pol}/rsd_per_feature.tsv")
        s = read_table(f"01_qc/{pol}/rsd_summary.tsv").set_index("mode")
        for mode, flag, col in [("peak_picked_only", "eligible_peak_picked",
                                 "rsd_peak_picked_pct"),
                                ("including_gap_filled", "eligible_with_gapfill",
                                 "rsd_with_gapfill_pct")]:
            v = np.sort(rsd.loc[rsd[flag].astype(bool), col].dropna().to_numpy())
            n = int(s.loc[mode, "n_eligible"])
            if len(v) != n:
                raise ValueError(f"{pol} {mode}: {len(v)} RSD values vs n_eligible {n}")
            vals.append(pd.DataFrame({"pol": pol, "mode": mode, "rsd_pct": v}))
            summ.append({"pol": pol, "mode": mode, "n_eligible": n,
                         "pct_le_30": s.loc[mode, "pct_of_eligible_le_30"],
                         "min_qc": int(s.loc[mode, "min_qc_values_required"]),
                         "n_qc": int(s.loc[mode, "n_pooled_qc"]),
                         "rsd_threshold_pct": thr["qc_feature_rsd_pct"]})
            number(f"qc_rsd_n_{mode}_{pol}", n, "qc_rsd")
            number(f"qc_rsd_pct_le_threshold_{mode}_{pol}",
                   float(s.loc[mode, "pct_of_eligible_le_30"]), "qc_rsd")
        number(f"qc_rsd_threshold_pct_{pol}", thr["qc_feature_rsd_pct"], "qc_rsd")
    table("qc_rsd_values", pd.concat(vals))
    table("qc_rsd_summary", pd.DataFrame(summ))
    number("qc_gapfill_feature_ratio_pos",
           NUMBERS["qc_rsd_n_including_gap_filled_pos"] / NUMBERS["qc_rsd_n_peak_picked_only_pos"],
           "qc_rsd")

    set_panel("qc_is")
    pts, info = [], []
    within, rt_ok = [], []
    for pol in POLS:
        thr = read_yaml(f"01_qc/{pol}/provenance.yaml")["params"]["thresholds_from_run"]
        ipi = read_table(f"01_qc/{pol}/is_per_injection.tsv")
        byb = read_table(f"01_qc/{pol}/is_by_block.tsv")
        ipi = ipi[ipi.detected.astype(bool)].merge(
            byb[["compound", "block", "median_rt"]], on=["compound", "block"], how="left")
        ipi["rt_minus_block_median"] = ipi.rt - ipi.median_rt
        wide = byb.pivot(index="compound", columns="block", values="median_rt").dropna()
        shift = (wide[2] - wide[1]).median()
        pts.append(ipi[["injection_order", "block", "compound", "rt_minus_block_median",
                        "ppm_error"]].assign(pol=pol))
        info.append({"pol": pol, "rt_tol_min": thr["rt_tolerance_min"],
                     "ppm_tol": thr["mass_accuracy_ppm"], "block_shift_min": shift,
                     "n_is_shift": len(wide), "n_is_detected": ipi.compound.nunique()})
        within.append(((ipi.ppm_error.abs() <= thr["mass_accuracy_ppm"]).sum(), len(ipi)))
        rt_ok.append(((ipi.rt_minus_block_median.abs() <= thr["rt_tolerance_min"]).sum(),
                      len(ipi)))
        number(f"qc_is_block_shift_min_{pol}", shift, "qc_is")
        number(f"qc_is_detected_{pol}", ipi.compound.nunique(), "qc_is")
        number(f"qc_rt_tolerance_min_{pol}", thr["rt_tolerance_min"], "qc_is")
        number(f"qc_mass_tolerance_ppm_{pol}", thr["mass_accuracy_ppm"], "qc_is")
    table("qc_is_points", pd.concat(pts))
    table("qc_is_info", pd.DataFrame(info))
    number("qc_is_pct_rt_within", sum(a for a, _ in rt_ok) / sum(b for _, b in rt_ok) * 100,
           "qc_is")
    number("qc_is_pct_ppm_within", sum(a for a, _ in within) / sum(b for _, b in within) * 100,
           "qc_is")

    # S4 (gated): internal standards against the neat calibration
    set_panel("qc_is_vs_cal")
    # m_found = extract concentration x reconstitution volume, so the exported amounts and
    # R_app scale with the volume; the export assumed recon_ul_ASSUMED, gates.yaml holds the
    # confirmed volume (Rahul, 2026-11-25: 200 uL).
    vol = gates().get("reconstitution_volume_ul")
    rows = []
    for pol in POLS:
        d = read_table(f"01_qc/{pol}/is_vs_calibration.tsv")
        assumed = float(read_yaml(f"01_qc/{pol}/provenance.yaml")["params"]["recon_ul_ASSUMED"])
        f = float(vol) / assumed if vol else 1.0
        for col in ("m_found_median_ng", "R_app_median_pct", "R_app_min_pct", "R_app_max_pct"):
            d[col] = d[col] * f
        rows.append(d.assign(pol=pol))
    d = pd.concat(rows)
    number("qc_recon_volume_ul", float(vol) if vol else assumed, "qc_is_vs_cal")
    number("figS4_shown", str(bool(gates().get("figS4_shown", True))), "qc_is_vs_cal")
    table("qc_is_vs_cal", d[d.n_samples_detected > 0])
    rec = d.dropna(subset=["R_app_median_pct"])
    table("qc_recovery", rec)
    number("qc_recovery_n_standards", len(rec), "qc_is_vs_cal")
    number("qc_recovery_median_min_pct", round(float(rec.R_app_median_pct.min()), 0), "qc_is_vs_cal")
    number("qc_recovery_median_max_pct", round(float(rec.R_app_median_pct.max()), 0), "qc_is_vs_cal")
    number("qc_recovery_n_above_100", int((rec.R_app_median_pct > 100).sum()), "qc_is_vs_cal")
    number("qc_is_n_within_cal", int((d.status == "within").sum()), "qc_is_vs_cal")
    number("qc_is_n_detected", int((d.n_samples_detected > 0).sum()), "qc_is_vs_cal")


# =========================================================================== annotation (S2)
SOURCE_ORDER = ["LipidSearch", "SIRIUS/CANOPUS", "MZmine lipid rules", "spectral library"]


def annotation():
    set_panel("ann_upset")
    sets_, inters = [], []
    for pol in POLS:
        s = read_table(f"02_annotation/{pol}/upset_sets.tsv")
        names = [x for x in SOURCE_ORDER if x in set(s.source)] + \
            [x for x in s.source if x not in SOURCE_ORDER]
        s = s.set_index("source").reindex(names).reset_index()
        sets_.append(s.assign(pol=pol, order=range(len(s))))
        i = read_table(f"02_annotation/{pol}/upset_intersections.tsv")
        inters.append(i.sort_values("n_features", ascending=False).assign(pol=pol))
    table("ann_upset_sets", pd.concat(sets_))
    table("ann_upset_intersections", pd.concat(inters))
    ip = inters[0]
    number("ann_single_source_pct_pos",
           ip[ip.n_sources == 1].n_features.sum() / ip.n_features.sum() * 100, "ann_upset")
    for pol, i in zip(POLS, inters):
        number(f"ann_features_with_source_{pol}", int(i.n_features.sum()), "ann_upset")

    set_panel("ann_confidence")
    rows = []
    for pol in POLS:
        c = read_table(f"02_annotation/{pol}/confidence.tsv").set_index("confidence")
        for lvl in ("Gold", "Silver", "Bronze"):
            rows.append({"pol": pol, "level": lvl, "n_features": int(c.loc[lvl, "n_features"]),
                         "n_all_features": int(c.n_features.sum())})
    conf = pd.DataFrame(rows)
    table("ann_confidence", conf)
    lv = conf[conf.pol == "pos"].set_index("level").n_features
    number("ann_top_level_pos", lv.idxmax(), "ann_confidence")
    number("ann_top_level_share_pos", lv.max() / lv.sum() * 100, "ann_confidence")
    number("ann_all_features_pos", int(conf[conf.pol == "pos"].n_all_features.iloc[0]),
           "ann_confidence")

    set_panel("ann_classes")
    rows = []
    for pol in POLS:
        c = read_table(f"02_annotation/{pol}/classes.tsv")
        rows.append(c.sort_values("n_features", ascending=False).head(12)
                    [["lipid_class", "n_features", "n_unique_lipids"]].assign(pol=pol))
        if pol == "pos":
            number("ann_top3_classes_pos", ", ".join(c.nlargest(3, "n_features").lipid_class),
                   "ann_classes")
            number("ann_adduct_redundancy_pct_pos",
                   (c.n_features.sum() - c.n_unique_lipids.sum()) / c.n_features.sum() * 100,
                   "ann_classes")
        else:
            number("ann_top_class_neg", c.nlargest(1, "n_features").lipid_class.iloc[0],
                   "ann_classes")
    table("ann_classes_top", pd.concat(rows))

    set_panel("ann_rt")
    rows = []
    for pol in POLS:
        d = read_table(f"02_annotation/{pol}/rt_model_PC.tsv")
        t = d[d.rt_check.isin(["ok", "outlier"])]
        rows.append(t[["rt_expected", "rt", "rt_check"]].assign(pol=pol))
        number(f"ann_rt_untested_{pol}", int((d.rt_check == "untested").sum()), "ann_rt")
        number(f"ann_rt_outliers_{pol}", int((d.rt_check == "outlier").sum()), "ann_rt")
        number(f"ann_rt_tested_{pol}", len(t), "ann_rt")
    table("ann_rt_pc", pd.concat(rows))


# =========================================================================== response (S3, S5)
ADDUCTS = {"M+H", "M+NH4", "M-H", "M+HCOO"}


def response():
    set_panel("resp_rie")
    d = read_table("03_response/response_by_class_adduct.tsv")
    if set(d.adduct) - ADDUCTS:
        raise ValueError(f"unmapped adducts: {set(d.adduct) - ADDUCTS}")
    order = d.groupby("class").logRIE_LPE.median().sort_values(ascending=False).index
    d = d.assign(class_rank=d["class"].map({c: i for i, c in enumerate(order)}))
    table("resp_rie", d[["class", "adduct", "RIE_LPE", "class_rank"]])
    number("resp_n_standards", len(d), "resp_rie")
    number("resp_rie_orders", np.log10(d.RIE_LPE.max() / d.RIE_LPE.min()), "resp_rie")

    set_panel("resp_ie")
    p = read_table("03_response/ie_loco_predictions.tsv")
    st = read_table("03_response/ie_model_stats.tsv")
    best = st[st.best.astype(bool)].iloc[0]
    table("resp_ie_loco", p[["pred_best", "logRIE_LPE"]])
    number("resp_ie_model", best.model, "resp_ie")
    number("resp_ie_r2", float(best.r2), "resp_ie")
    number("resp_ie_rmse_log10", float(best.rmse_log10), "resp_ie")
    number("resp_ie_fold_error", float(best.fold_error), "resp_ie")
    number("resp_ie_n_standards", int(best.n_standards), "resp_ie")
    number("resp_ie_n_classes", int(best.n_classes), "resp_ie")

    set_panel("resp_nist")
    rows = []
    for pol in POLS:
        n = read_table(f"05_nist/{pol}/nist_comparison.tsv")
        n = n.assign(pol=pol, route_group=np.where(n.routes == "Q1", "Q1", "Q3"))
        rows.append(n[["pol", "route_group", "medm", "estimate_nmol_per_ml", "ratio",
                       "within_2fold"]])
        q1 = n[n.routes == "Q1"]
        number(f"nist_q1_n_{pol}", len(q1), "resp_nist")
        if len(q1):
            number(f"nist_q1_median_ratio_{pol}", float(np.median(q1.ratio)), "resp_nist")
            number(f"nist_q1_pct_within_2fold_{pol}",
                   q1.within_2fold.astype(bool).mean() * 100, "resp_nist")
        number(f"nist_q3_n_{pol}", int((n.routes != "Q1").sum()), "resp_nist")
        q3 = n[n.routes != "Q1"]
        if len(q3):
            number(f"nist_q3_median_ratio_{pol}", float(np.median(q3.ratio)), "resp_nist")
            number(f"nist_q3_pct_within_2fold_{pol}",
                   q3.within_2fold.astype(bool).mean() * 100, "resp_nist")
    table("resp_nist", pd.concat(rows))

    # S5 (gated): class amounts per gram, common lipid set
    set_panel("quant_class_amounts")
    summ = []
    for pol in POLS:
        d = read_table(f"04_quant/{pol}/class_by_sample.tsv")
        g = d.groupby("lipid_class").ug_per_g.agg(["mean", "std", "count"]).reset_index()
        summ.append(g.assign(pol=pol))
    table("quant_class_amounts", pd.concat(summ))


# =========================================================================== Atlas (Fig 3, GA)
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


def atlas():
    set_panel("atlas")
    cov = read_table("06_atlas_match/pos/coverage.tsv").iloc[0]
    conf = read_table("02_annotation/pos/confidence.tsv")
    prov = read_yaml("06_atlas_match/pos/provenance.yaml")["params"]
    det = read_table("06_atlas_match/pos/phylum_detection.tsv")
    m = _matched()
    n_ms2, n_match = int(cov.soil_features_with_ms2), int(cov.soil_features_matched)
    if n_match != len(m):
        raise ValueError(f"coverage {n_match} vs best_match rows {len(m)}")
    n_named = int(m.named.sum())
    number("atlas_n_features_all", int(conf.n_features.sum()), "atlas")
    number("atlas_n_ms2", n_ms2, "atlas")
    number("atlas_n_matched", n_match, "atlas")
    number("atlas_n_matched_named", n_named, "atlas")
    number("atlas_n_matched_unnamed", n_match - n_named, "atlas")
    number("atlas_pct_matched", float(cov.pct_matched), "atlas")
    number("atlas_pct_unnamed_of_matched", (n_match - n_named) / n_match * 100, "atlas")
    number("atlas_pct_single_group", (m.n_groups == 1).mean() * 100, "atlas")
    number("atlas_median_ref_phyla", float(m.n_ref_phyla.median()), "atlas")
    number("atlas_match_min_cosine", prov["min_cosine"], "atlas")
    number("atlas_match_min_peaks", prov["min_matched_peaks"], "atlas")
    number("atlas_match_ppm", prov["ppm"], "atlas")
    grp = atlas_groups().set_index("phylum").group
    units = det.groupby("phylum").n_units_phylum.first()
    number("atlas_n_phyla_matching", len(units), "atlas")
    number("atlas_n_ref_samples_matching", int(units.sum()), "atlas")

    table("atlas_coverage", pd.DataFrame([
        {"segment": "matched, no name", "n": n_match - n_named},
        {"segment": "matched and named", "n": n_named},
        {"segment": "no Atlas match", "n": n_ms2 - n_match}]))

    bins = np.arange(1, len(GROUP_ORDER) + 1)
    table("atlas_breadth", pd.DataFrame({
        "n_groups": bins,
        "unnamed": m[~m.named].n_groups.value_counts().reindex(bins, fill_value=0).to_numpy(),
        "named": m[m.named].n_groups.value_counts().reindex(bins, fill_value=0).to_numpy()}))

    samples = units.groupby(grp).sum()
    occ = m.explode("groups").groups.value_counts().reindex(GROUP_ORDER, fill_value=0)
    only = m[m.n_groups == 1].groups.str[0].value_counts().reindex(GROUP_ORDER, fill_value=0)
    table("atlas_per_group", pd.DataFrame({
        "group": GROUP_ORDER, "ref_samples": [int(samples[g]) for g in GROUP_ORDER],
        "occurs_pct": [occ[g] / len(m) * 100 for g in GROUP_ORDER],
        "only_pct": [only[g] / len(m) * 100 for g in GROUP_ORDER]}))

    a = m[m.named]
    order = a.lipid_class.value_counts().head(8).index.tolist()
    e = a[a.lipid_class.isin(order)].explode("groups")
    share = (e.groupby(["lipid_class", "groups"]).soil_feature_id.nunique().unstack(fill_value=0)
             .reindex(index=order, columns=GROUP_ORDER, fill_value=0)
             .div(a.lipid_class.value_counts()[order], axis=0) * 100)
    cg = share.stack().rename("share_pct").reset_index().rename(columns={"groups": "group"})
    cg["n_features"] = cg.lipid_class.map(a.lipid_class.value_counts())
    cg["class_rank"] = cg.lipid_class.map({c: i for i, c in enumerate(order)})
    table("atlas_class_group", cg)

    # graphical abstract: release reference samples per organism group (analysed phyla)
    set_panel("atlas_release")
    g = atlas_groups()
    g = g[g.analysed.astype(str) == "True"]
    per = g.groupby("group").samples_pos.sum().reindex(GROUP_ORDER)
    table("atlas_release_groups", pd.DataFrame({"group": GROUP_ORDER,
                                                "ref_samples": per.astype(int).to_numpy()}))
    number("atlas_release_ref_samples", int(per.sum()), "atlas_release")
    number("atlas_release_phyla", len(g), "atlas_release")


# =========================================================================== interpretation (Fig 4)
def _shared_fragments(soil, atlas_, n=4):
    out = []
    a_mz = atlas_.mz.to_numpy()
    for _, r in soil.iterrows():
        j = np.argmin(np.abs(a_mz - r.mz))
        if abs(a_mz[j] - r.mz) <= FRAG_TOL_DA:
            out.append((min(r.rel_intensity, atlas_.rel_intensity.iloc[j]), r.mz,
                        r.rel_intensity))
    return sorted(out, reverse=True)[:n]


def interpretation():
    set_panel("interp")
    exs = read_table("06_atlas_match/examples/examples.tsv")
    table("interp_examples", exs)
    grp = atlas_groups().set_index("phylum").group
    det = read_table("06_atlas_match/examples/examples_phylum_detection.tsv")
    for role in ("broad", "restricted"):
        ex = exs[exs.role == role].iloc[0]
        m = read_table(f"06_atlas_match/examples/mirror_{ex.soil_feature_id}.tsv")
        soil, atl = m[m.side == "soil"], m[m.side == "atlas"]
        xmax = ex.soil_mz + 25
        placed, labels = [], []
        for _, mz, h in _shared_fragments(soil, atl):
            if any(abs(mz - p) < 45 for p in placed):  # keep labels from overlapping
                continue
            placed.append(mz)
            labels.append({"mz": mz, "height": h})
        table(f"interp_mirror_{role}", m[["mz", "rel_intensity", "side"]])
        table(f"interp_labels_{role}", pd.DataFrame(labels, columns=["mz", "height"]))
        number(f"interp_xmax_{role}", xmax, "interp")
        number(f"interp_hidden_peaks_{role}", int((m.mz > xmax).sum()), "interp")
        gs = sorted({grp[p] for p in det[det.feature_id == ex.ref_feature_id].phylum})
        number(f"interp_groups_{role}", ";".join(gs), "interp")
    units_all = read_table("06_atlas_match/pos/phylum_detection.tsv") \
        .groupby("phylum").n_units_phylum.first()
    cols = [(g, p) for g in GROUP_ORDER for p in grp[grp == g].index if p in units_all.index]
    rows = []
    for role in ("broad", "restricted"):
        ex = exs[exs.role == role].iloc[0]
        d = det[det.feature_id == ex.ref_feature_id].set_index("phylum")
        for j, (g, p) in enumerate(cols):
            hit = p in d.index
            rows.append({"role": role, "col": j, "group": g, "phylum": p, "detected": hit,
                         "n_detected": int(d.loc[p, "n_detected"]) if hit else 0,
                         "n_units": int(d.loc[p, "n_units_phylum"]) if hit
                         else int(units_all[p]),
                         "detection_rate": float(d.loc[p, "detection_rate"]) if hit else 0.0})
    table("interp_occurrence", pd.DataFrame(rows))


# =========================================================================== application (Fig 5, S5)
def application():
    set_panel("app")
    s = read_table("00_inputs/climgrass_samples.csv")
    table("app_design", s[["plot", "treatment", "warming_co2", "drought", "used_pos", "used_neg"]])
    number("app_soils_pos", int((s.used_pos == "yes").sum()), "app")
    number("app_soils_neg", int((s.used_neg == "yes").sum()), "app")
    resp = {p: read_table(f"08_application/{p}/responsive.tsv") for p in POLS}
    rows = pd.concat([resp[p].assign(pol=p) for p in POLS]).sort_values(
        "log2fc_drought_vs_ambient")
    table("app_responsive", rows[["pol", "feature_id", "lipid_name",
                                  "log2fc_drought_vs_ambient"]])
    for p in POLS:
        ts = read_table(f"08_application/{p}/test_summary.tsv").iloc[0]
        number(f"app_n_lipids_tested_{p}", int(ts.n_lipids_tested), "app")
        number(f"app_n_q05_climate_{p}", int(ts["n_q_lt_0.05_climate"]), "app")
        number(f"app_n_q05_interaction_{p}", int(ts["n_q_lt_0.05_interaction"]), "app")
    for p in POLS:
        number(f"app_n_responsive_{p}", len(resp[p]), "app")
    occ = read_table("06_atlas_match/examples/responsive_atlas_occurrence.tsv")
    number("app_n_responsive_pos_matched", int(occ.n_ref_phyla.notna().sum()), "app")
    table("app_occurrence", occ.sort_values("n_ref_phyla", na_position="first")
          [["feature_id", "lipid_name", "n_ref_phyla", "cosine"]])
    trace = resp["pos"].merge(occ[["feature_id", "cosine", "n_ref_phyla"]], on="feature_id",
                              how="left")
    # organism groups of the reference phyla of each responsive POS lipid (Fig 5D)
    grp = atlas_groups()
    grp = grp[grp.analysed.astype(str) == "True"]
    n_an = grp.groupby("group").phylum.size()
    gmap = grp.set_index("phylum").group
    rows_g = []
    for _, r in occ.iterrows():
        ref = set(str(r.ref_phyla).split(";")) if pd.notna(r.ref_phyla) else set()
        for g in GROUP_ORDER:
            rows_g.append({"lipid_name": r.lipid_name, "group": g,
                           "n_phyla_with_ref": sum(gmap.get(x) == g for x in ref),
                           "n_phyla_analysed": int(n_an.get(g, 0))})
    table("app_trace_groups", pd.DataFrame(rows_g))
    lipid = gates().get("fig7E_lipid")  # chosen by Rahul in gates.yaml; drawn as Fig 5D
    number("app_trace_lipid", lipid or "none", "app")
    table("app_trace", trace[["lipid_name", "confidence", "quant_route", "q_drought", "cosine",
                              "n_ref_phyla", "mean_ng_g_ambient", "mean_ng_g_future",
                              "mean_ng_g_ambient_drought", "mean_ng_g_future_drought"]])
    set_panel("app_profiles")
    table("app_class_profiles", read_table("08_application/pos/class_profiles.tsv"))
    read_table("08_application/neg/class_profiles.tsv")


# =========================================================================== run
def run() -> None:
    reset_reads()
    TABLES.clear(), NUMBERS.clear(), SOURCES.clear()
    for section in (release, qc, annotation, response, atlas, interpretation, application):
        section()
    DERIVED.mkdir(exist_ok=True)
    for old in DERIVED.glob("*.tsv"):
        old.unlink()
    for name, df in TABLES.items():
        df.to_csv(DERIVED / f"{name}.tsv", sep="\t", index=False)
    pd.DataFrame({"key": list(NUMBERS), "value": [NUMBERS[k] for k in NUMBERS]}).to_csv(
        DERIVED / "numbers.tsv", sep="\t", index=False)
    (DERIVED / "sources.yaml").write_text(yaml.safe_dump(SOURCES, sort_keys=True))


def _same(a: str, b: str) -> bool:
    try:
        return abs(float(a) - float(b)) <= 1e-9 * max(1.0, abs(float(b)))
    except ValueError:
        return a == b


def check() -> list[str]:
    """Compare numbers.tsv with expected_numbers.csv; return the differences."""
    got = pd.read_csv(DERIVED / "numbers.tsv", sep="\t", dtype=str).set_index("key").value
    exp = pd.read_csv(EXPECTED, dtype=str).set_index("key").value
    diffs = [f"missing: {k}" for k in exp.index if k not in got.index]
    diffs += [f"new, not expected: {k}" for k in got.index if k not in exp.index]
    diffs += [f"{k}: expected {exp[k]}, got {got[k]}" for k in exp.index
              if k in got.index and not _same(got[k], exp[k])]
    return diffs


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--check", action="store_true", help="compare with expected_numbers.csv")
    ap.add_argument("--record", action="store_true",
                    help="write expected_numbers.csv from this run (only after review)")
    a = ap.parse_args()
    run()
    print(f"wrote {len(TABLES)} tables and {len(NUMBERS)} numbers to {DERIVED}")
    if a.record:
        pd.read_csv(DERIVED / "numbers.tsv", sep="\t", dtype=str).to_csv(EXPECTED, index=False)
        print(f"recorded {EXPECTED}")
    if a.check:
        diffs = check()
        if diffs:
            print("CHECK FAILED\n  " + "\n  ".join(diffs))
            return 1
        print(f"check passed: {len(NUMBERS)} numbers match expected_numbers.csv")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
