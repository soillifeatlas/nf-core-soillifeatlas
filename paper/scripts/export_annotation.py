#!/usr/bin/env python3
"""Issue #2: export the merged ClimGrass annotation table for Fig. 4 and section 3.4.

Writes per polarity:
  annotations.tsv        per feature: sources, confidence, class, composition_check, rt residual, flags
  upset_sets.tsv         set size per annotation source
  upset_intersections.tsv  exclusive intersection counts (each feature counted once)
  confidence.tsv         Gold / Silver / Bronze / unidentified (denominator = all MZmine features)
  classes.tsv            per class: features and unique lipids (adduct redundancy collapsed)
  rt_model_<class>.tsv   rt, carbons, double bonds, expected rt and flag for one class
Unique lipid rule: Gold/Silver rows with the same lipid_name whose RTs lie within
--collapse-rt min of each other are one lipid (different adducts of one molecule).
Bronze (class-only) rows cannot be collapsed and count one each.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from soillifeatlas import annotate as an  # noqa: E402
import provenance  # noqa: E402


def collapse_unique(ann: pd.DataFrame, rt_window: float) -> pd.Series:
    """Return a lipid id per row; rows with the same name within rt_window share it."""
    ids = pd.Series(index=ann.index, dtype=object)
    named = ann["confidence"].isin(["Gold", "Silver"])
    for name, g in ann[named].sort_values("rt").groupby("lipid_name"):
        k, last = 0, None
        for idx, rt in g["rt"].items():
            if last is not None and rt - last > rt_window:
                k += 1
            ids[idx] = f"{name}#{k}"
            last = rt
    ids[~named] = "feature:" + ann.loc[~named, "feature_id"].astype(str)
    return ids


def upset(ann: pd.DataFrame):
    sets = ann["sources"].str.split("; ")
    sizes = sets.explode().value_counts().rename_axis("source").reset_index(name="n_features")
    inter = ann["sources"].value_counts().rename_axis("intersection").reset_index(name="n_features")
    inter["n_sources"] = inter["intersection"].str.count(";") + 1
    return sizes, inter


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--results", required=True)
    ap.add_argument("--prefix", required=True)
    ap.add_argument("--run-provenance", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--collapse-rt", type=float, default=0.1)
    ap.add_argument("--rt-class", default=None, help="class for the RT-model panel (default: most Gold anchors)")
    a = ap.parse_args(argv)
    res, out = Path(a.results), Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    ann_p = res / "annotate" / f"{a.prefix}.annotations.tsv"
    full_p = res / "mzmine" / f"{a.prefix}.features_full.csv"
    ann = pd.read_csv(ann_p, sep="\t")
    n_features = sum(1 for _ in open(full_p)) - 1

    ann["rt_residual_min"] = (ann["rt"] - ann["rt_expected"]).round(3)
    ann["lipid_id"] = collapse_unique(ann, a.collapse_rt)
    flags = []
    for _, r in ann.iterrows():
        f = []
        if r["composition_check"] == "implausible":
            f.append("composition_implausible")
        if r["rt_check"] == "outlier":
            f.append("rt_outlier")
        if not r["class_agreement"]:
            f.append("class_disagreement")
        flags.append(";".join(f))
    ann["flags"] = flags
    ann.to_csv(out / "annotations.tsv", sep="\t", index=False)

    sizes, inter = upset(ann)
    sizes.to_csv(out / "upset_sets.tsv", sep="\t", index=False)
    inter.to_csv(out / "upset_intersections.tsv", sep="\t", index=False)

    conf = ann["confidence"].value_counts().reindex(["Gold", "Silver", "Bronze"]).fillna(0).astype(int)
    conf["unidentified"] = n_features - len(ann)
    c = conf.rename_axis("confidence").reset_index(name="n_features")
    c["pct_of_all_features"] = (c["n_features"] / n_features * 100).round(2)
    c["n_unique_lipids"] = [ann.loc[ann["confidence"] == k, "lipid_id"].nunique() if k != "unidentified" else np.nan
                            for k in c["confidence"]]
    c.to_csv(out / "confidence.tsv", sep="\t", index=False)

    cls = ann.groupby("lipid_class").agg(
        n_features=("feature_id", "size"),
        n_unique_lipids=("lipid_id", "nunique"),
        n_gold=("confidence", lambda s: int((s == "Gold").sum())),
        n_silver=("confidence", lambda s: int((s == "Silver").sum())),
        n_bronze=("confidence", lambda s: int((s == "Bronze").sum())),
        n_rt_outlier=("rt_check", lambda s: int((s == "outlier").sum())),
        n_composition_implausible=("composition_check", lambda s: int((s == "implausible").sum())),
    ).sort_values("n_features", ascending=False)
    cls.to_csv(out / "classes.tsv", sep="\t")

    tested = ann[ann["rt_check"].isin(["ok", "outlier"])]
    rt_class = a.rt_class or (tested[tested["confidence"] == "Gold"]["lipid_class"].value_counts().idxmax()
                              if len(tested) else None)
    if rt_class:
        g = ann[ann["lipid_class"] == rt_class].copy()
        cd = g["lipid_name"].map(an._carbons_dbs)
        g["carbons"], g["double_bonds"] = [x[0] for x in cd], [x[1] for x in cd]
        g["rt_group"] = [an._rt_group(n, k) for n, k in zip(g["lipid_name"], g["lipid_class"])]
        g["is_anchor"] = (g["confidence"] == "Gold") & g["carbons"].notna() & (g["composition_check"] != "implausible")
        g[["feature_id", "lipid_name", "confidence", "rt_group", "carbons", "double_bonds", "rt", "rt_expected",
           "rt_residual_min", "rt_check", "is_anchor"]].to_csv(out / f"rt_model_{rt_class}.tsv", sep="\t", index=False)

    summary = {
        "n_mzmine_features": n_features, "n_annotated_features": len(ann),
        "n_unique_lipids_gold_silver": int(ann.loc[ann["confidence"].isin(["Gold", "Silver"]), "lipid_id"].nunique()),
        "n_rt_outlier": int((ann["rt_check"] == "outlier").sum()),
        "n_rt_tested": int(ann["rt_check"].isin(["ok", "outlier"]).sum()),
        "n_composition_implausible": int((ann["composition_check"] == "implausible").sum()),
        "rt_model_class": rt_class, "quant_route_in_annotate": ann["quant_route"].value_counts().to_dict(),
    }
    pd.Series(summary).to_csv(out / "summary.tsv", sep="\t", header=False)
    provenance.write(out, {"annotations": ann_p, "features_full": full_p},
                     {"collapse_rt_min": a.collapse_rt, "rt_class": rt_class},
                     pipeline_run=provenance.read_run_provenance(a.run_provenance))


if __name__ == "__main__":
    main()
