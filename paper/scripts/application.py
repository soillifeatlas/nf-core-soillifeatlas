#!/usr/bin/env python3
"""Issue #8: ClimGrass application (Fig. 7B-C, section 3.6).

From the stage-5 output of one polarity:
  class_profiles.tsv     ug per g dry soil per class and treatment: mean, SD, n (plots)
  lipid_tests.tsv        per unique lipid (Q1-Q3): two-way ANOVA on log10(ng/g),
                         factors drought (no/yes) x future climate (no/yes), type II SS;
                         BH-FDR across lipids per effect; fold change vs ambient
  responsive.tsv         lipids with q < --fdr for drought, climate or the interaction
Only lipids quantified in >= 2 plots of every treatment are tested.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import provenance  # noqa: E402

TREAT = {"ambient": (0, 0), "future": (0, 1), "ambient_drought": (1, 0), "future_drought": (1, 1)}


def bh(p: pd.Series) -> pd.Series:
    p = p.astype(float)
    ok = p.notna()
    s = p[ok].sort_values(ascending=False)
    rank = np.arange(len(s), 0, -1)
    q = pd.Series(np.minimum.accumulate(s.to_numpy() * len(s) / rank), index=s.index).clip(upper=1)
    return q.reindex(p.index)


def anova2(y: np.ndarray, d: np.ndarray, c: np.ndarray) -> dict:
    """Type II two-way ANOVA with interaction for a 2x2 design (numpy least squares)."""
    from scipy import stats

    def rss(cols):
        X = np.column_stack([np.ones(len(y))] + cols)
        beta, *_ = np.linalg.lstsq(X, y, rcond=None)
        return float(((y - X @ beta) ** 2).sum()), X.shape[1]

    full, kf = rss([d, c, d * c])
    df_res = len(y) - kf
    if df_res <= 0:
        return {}
    mse = full / df_res
    add, _ = rss([d, c])
    only_c, _ = rss([c])
    only_d, _ = rss([d])
    out = {}
    for name, ss in (("drought", only_c - add), ("climate", only_d - add), ("interaction", add - full)):
        F = ss / mse if mse > 0 else np.nan
        out[f"F_{name}"] = F
        out[f"p_{name}"] = float(stats.f.sf(F, 1, df_res)) if np.isfinite(F) else np.nan
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--quant", required=True, help="paper/04_quant/<pol>")
    ap.add_argument("--out", required=True)
    ap.add_argument("--fdr", type=float, default=0.05)
    a = ap.parse_args(argv)
    q, out = Path(a.quant), Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    est_p, cls_p, per_p = q / "estimates_long.tsv", q / "class_by_sample.tsv", q / "per_lipid.tsv"

    cls = pd.read_csv(cls_p, sep="\t")
    prof = cls.groupby(["lipid_class", "treatment"])["ug_per_g"].agg(["mean", "std", "count"]).reset_index()
    prof.rename(columns={"count": "n_plots"}).to_csv(out / "class_profiles.tsv", sep="\t", index=False)
    tot = cls.groupby(["plot", "treatment"])["ug_per_g"].sum().reset_index()
    tot.groupby("treatment")["ug_per_g"].agg(["mean", "std", "count"]).to_csv(out / "total_by_treatment.tsv", sep="\t")

    est = pd.read_csv(est_p, sep="\t")
    per = pd.read_csv(per_p, sep="\t")
    est = est[est["quant_route"].isin(["Q1", "Q2", "Q3"])].merge(per[["feature_id", "lipid_id", "mz", "rt"]],
                                                                    on="feature_id")
    # one feature per unique lipid: largest median area
    keep = est.groupby(["lipid_id", "feature_id"])["area"].median().reset_index() \
              .sort_values("area", ascending=False).drop_duplicates("lipid_id")["feature_id"]
    est = est[est["feature_id"].isin(keep) & est["ng_per_g"].gt(0)]
    rows = []
    for fid, g in est.groupby("feature_id"):
        n_cell = g.groupby("treatment")["ng_per_g"].count().reindex(list(TREAT)).fillna(0)
        if (n_cell < 2).any():
            continue
        d = g["treatment"].map(lambda t: TREAT[t][0]).to_numpy(float)
        c = g["treatment"].map(lambda t: TREAT[t][1]).to_numpy(float)
        res = anova2(np.log10(g["ng_per_g"].to_numpy(float)), d, c)
        if not res:
            continue
        m = g.groupby("treatment")["ng_per_g"].mean()
        r0 = g.iloc[0]
        rows.append({"feature_id": fid, "lipid_name": r0["lipid_name"], "lipid_class": r0["lipid_class"],
                     "confidence": r0["confidence"], "quant_route": r0["quant_route"], "mz": r0["mz"], "rt": r0["rt"],
                     **{f"mean_ng_g_{t}": m.get(t, np.nan) for t in TREAT},
                     "log2fc_drought_vs_ambient": np.log2(m.get("ambient_drought") / m.get("ambient")),
                     "log2fc_future_vs_ambient": np.log2(m.get("future") / m.get("ambient")),
                     **res})
    t = pd.DataFrame(rows)
    for eff in ("drought", "climate", "interaction"):
        t[f"q_{eff}"] = bh(t[f"p_{eff}"])
    t.to_csv(out / "lipid_tests.tsv", sep="\t", index=False)
    resp = t[(t[["q_drought", "q_climate", "q_interaction"]] < a.fdr).any(axis=1)]
    resp.sort_values("q_drought").to_csv(out / "responsive.tsv", sep="\t", index=False)
    pd.DataFrame([{"n_lipids_tested": len(t), **{f"n_q_lt_{a.fdr}_{e}": int((t[f"q_{e}"] < a.fdr).sum())
                                                  for e in ("drought", "climate", "interaction")}}]).to_csv(
        out / "test_summary.tsv", sep="\t", index=False)
    provenance.write(out, {"estimates_long": est_p, "class_by_sample": cls_p, "per_lipid": per_p},
                     {"test": "two-way ANOVA (type II) on log10 ng/g, drought x future climate, n = 3 plots per cell",
                      "correction": "Benjamini-Hochberg per effect", "fdr": a.fdr,
                      "inclusion": ">= 2 quantified plots in every treatment; one feature per unique lipid"})


if __name__ == "__main__":
    main()
