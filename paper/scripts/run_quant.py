#!/usr/bin/env python3
"""Issue #4: run stage 5 (soillifeatlas.quantify) on one ClimGrass polarity.

Writes to --out:
  estimates_long.tsv    feature x soil: ng per g dry soil, route, standard, f, f source, range flag
  per_lipid.tsv         one row per feature: route, standard, f, f source, flags, median ng/g
  route_counts.tsv      features and unique lipids per route
  class_by_sample.tsv   ug per g dry soil per class and soil (Q1-Q3, one feature per unique lipid)
  is_range.tsv          standard area in each soil against its neat calibration range
Areas: peak-picked values only (gap-filled values are not quantified).
"""
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from soillifeatlas import process_qc as pq  # noqa: E402
from soillifeatlas import quantify as qt  # noqa: E402
import provenance  # noqa: E402
from export_annotation import collapse_unique  # noqa: E402


def plot_of(file_id: str):
    m = re.match(r"OE12-3-(\d+)(POS|NEG)$", file_id)
    return int(m.group(1)) if m else None


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--results", required=True)
    ap.add_argument("--prefix", required=True)
    ap.add_argument("--polarity", choices=["pos", "neg"], required=True)
    ap.add_argument("--samples", required=True, help="paper/00_inputs/climgrass_samples.csv")
    ap.add_argument("--rie", required=True, help="Samrat et al. 2025 Table S10 (rie_table_s10.csv)")
    ap.add_argument("--predicted", default=None, help="optional lipid_name,RIE_pred table (Q2)")
    ap.add_argument("--run-provenance", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--collapse-rt", type=float, default=0.1)
    a = ap.parse_args(argv)
    res, out = Path(a.results), Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    p = {"features_full": res / "mzmine" / f"{a.prefix}.features_full.csv",
         "annotations": res / "annotate" / f"{a.prefix}.annotations.tsv",
         "internal_standards": res / "qc" / f"{a.prefix}.qc_internal_standards.tsv",
         "is_per_injection": res / "qc" / f"{a.prefix}.qc_is_per_injection.tsv",
         "samples": Path(a.samples), "rie": Path(a.rie)}
    if a.predicted:
        p["predicted"] = Path(a.predicted)

    full = pd.read_csv(p["features_full"], low_memory=False)
    area = pq.per_file_matrix(full, "area")
    area = area.where(pq.detected_matrix(full, area))           # peak-picked only
    ann = pd.read_csv(p["annotations"], sep="\t")
    ann["lipid_id"] = collapse_unique(ann, a.collapse_rt)

    ist = pd.read_csv(p["internal_standards"], sep="\t")
    ist = ist[ist["matched"].astype(str).str.lower() == "true"]
    standards = ist.rename(columns={"median_rt": "rt"})[["compound", "lipid_class", "rt", "adduct"]]
    ipi = pd.read_csv(p["is_per_injection"], sep="\t")
    det = ipi[ipi["detected"] & (ipi["area"] > 0)]
    is_area = det[det["role"] == "sample"].pivot_table(index="compound", columns="file_id", values="area")
    cal = det[det["role"] == "standard_cal"].groupby("compound")["area"].agg(cal_min_area="min",
                                                                               cal_max_area="max").reset_index()

    meta = pd.read_csv(p["samples"])
    use = "used_pos" if a.polarity == "pos" else "used_neg"
    meta = meta[meta[use] == "yes"]
    soil_files = [f for f in area.columns if plot_of(f) in set(meta["plot"])]
    samples = pd.DataFrame({"file_id": soil_files, "plot": [plot_of(f) for f in soil_files]}).merge(
        meta[["plot", "treatment", "dry_mass_g", "spike_ng_per_is"]], on="plot")
    samples = samples.rename(columns={"spike_ng_per_is": "spike_ng"})

    rie = pd.read_csv(p["rie"])
    pred = pd.read_csv(p["predicted"]) if a.predicted else None
    routes = qt.assign_routes(ann, standards, pred)
    est = qt.estimate(routes, area, is_area, samples, standards, rie, cal, a.polarity)
    est = est.merge(samples[["file_id", "plot", "treatment"]], on="file_id")
    est.to_csv(out / "estimates_long.tsv", sep="\t", index=False)

    per = routes.merge(ann[["feature_id", "mz", "rt", "lipid_id"]], on="feature_id")
    if len(est):
        agg = est.groupby("feature_id").agg(
            f=("f", "first"), f_source=("f_source", "first"), f_default=("f_default", "first"),
            n_soils_quantified=("ng_per_g", lambda s: int(s.notna().sum())),
            median_ng_per_g=("ng_per_g", "median"),
            is_range_flags=("is_range_flag", lambda s: ";".join(sorted(set(s))))).reset_index()
        per = per.merge(agg, on="feature_id", how="left")
    per.to_csv(out / "per_lipid.tsv", sep="\t", index=False)

    n_feat_total = len(area)
    rc = per.groupby("quant_route").agg(n_features=("feature_id", "size"),
                                        n_unique_lipids=("lipid_id", "nunique")).reindex(qt.ROUTE_ORDER).fillna(0)
    rc.loc["unannotated"] = [n_feat_total - len(per), np.nan]
    rc["pct_of_annotated_features"] = (rc["n_features"] / len(per) * 100).round(1)
    rc.loc["unannotated", "pct_of_annotated_features"] = np.nan
    rc.to_csv(out / "route_counts.tsv", sep="\t")

    # class sums: one feature per unique lipid (the one with the largest median area), Q1-Q3 only
    q = est[est["quant_route"].isin(["Q1", "Q2", "Q3"])].merge(ann[["feature_id", "lipid_id"]], on="feature_id")
    if len(q):
        best = q.groupby(["lipid_id", "feature_id"])["area"].median().reset_index() \
                .sort_values("area", ascending=False).drop_duplicates("lipid_id")["feature_id"]
        q = q[q["feature_id"].isin(best)]
        cls = q.groupby(["lipid_class", "plot", "treatment"])["ng_per_g"].sum(min_count=1).reset_index()
        cls["ug_per_g"] = cls["ng_per_g"] / 1000
        cls.to_csv(out / "class_by_sample.tsv", sep="\t", index=False)

    rng = is_area.reindex(columns=soil_files).stack().rename("is_area").reset_index()
    rng = rng.merge(cal, on="compound", how="left")
    rng["flag"] = [qt.range_flag(r.is_area, r.cal_min_area, r.cal_max_area) for r in rng.itertuples()]
    rng.to_csv(out / "is_range.tsv", sep="\t", index=False)

    provenance.write(out, p, {"polarity": a.polarity, "areas": "peak-picked only",
                              "collapse_rt_min": a.collapse_rt, "q2_predictions": bool(a.predicted),
                              "routes": "option A (manuscript 2.7), decided 2026-11-25"},
                     notes="Calibration range = min/max neat ESTD-1 areas (sequence OE12-1, 27-28 Mar 2023).",
                     pipeline_run=provenance.read_run_provenance(a.run_provenance))


if __name__ == "__main__":
    main()
