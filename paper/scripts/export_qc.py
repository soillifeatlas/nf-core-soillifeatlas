#!/usr/bin/env python3
"""Issue #1: export the ClimGrass QC tables for Fig. 3 and section 3.3.

Reads one polarity's pipeline outputs (results/mzmine, results/qc) and writes:
  injections.tsv          injection order, role, block, start time
  rsd_per_feature.tsv     pooled-QC RSD peak-picked only vs including gap-filled
  rsd_summary.tsv         eligible-feature counts, denominators, medians, share <= 30 %
  is_per_injection.tsv    IS RT, RT deviation, mass error and area per injection, with block
  is_vs_calibration.tsv   soil IS amounts against the neat ESTD-1 calibration (in/below/above)
  blank_summary.tsv       sample features also high in blanks
  drift.tsv, checks.tsv   copies of the stage-3 tables
Thresholds come from the run's qc_summary.json, not from this script.
"""
from __future__ import annotations

import argparse
import json
import shutil
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from soillifeatlas import process_qc as pq  # noqa: E402
import provenance  # noqa: E402

# ESTD-1 (Notion "ESTD-Aliquot-1"): EquiSPLASH 100 ng/uL per component, 1:10 to level A1,
# then 2-fold serial dilution to A7. Nominal vial concentration in ng/uL.
ESTD1_NG_PER_UL = {f"A{i}": 10.0 / 2 ** (i - 1) for i in range(1, 8)}


def rsd_tables(area: pd.DataFrame, det: pd.DataFrame, qcs: list[str], min_frac: float, limit: float):
    need = max(2, int(np.ceil(min_frac * len(qcs))))
    qa = area[qcs]
    picked = qa.where(det[qcs] & (qa > 0))
    allv = qa.where(qa > 0)
    n_picked, n_all = picked.notna().sum(axis=1), allv.notna().sum(axis=1)

    def rsd(m):
        return m.std(axis=1, ddof=1) / m.mean(axis=1) * 100

    per = pd.DataFrame({
        "qc_peak_picked": n_picked, "qc_with_value": n_all,
        "eligible_peak_picked": n_picked >= need, "eligible_with_gapfill": n_all >= need,
        "rsd_peak_picked_pct": rsd(picked), "rsd_with_gapfill_pct": rsd(allv),
    })
    per.loc[~per["eligible_peak_picked"], "rsd_peak_picked_pct"] = np.nan
    per.loc[~per["eligible_with_gapfill"], "rsd_with_gapfill_pct"] = np.nan
    rows = []
    for mode, col, elig in (("peak_picked_only", "rsd_peak_picked_pct", "eligible_peak_picked"),
                            ("including_gap_filled", "rsd_with_gapfill_pct", "eligible_with_gapfill")):
        r = per.loc[per[elig], col]
        rows.append({"mode": mode, "n_features_total": len(per), "n_pooled_qc": len(qcs),
                     "min_qc_values_required": need, "n_eligible": int(per[elig].sum()),
                     "median_rsd_pct": round(float(r.median()), 2),
                     f"n_rsd_le_{limit:g}": int((r <= limit).sum()),
                     f"pct_of_eligible_le_{limit:g}": round(float((r <= limit).mean() * 100), 1)})
    # same features in both modes: the effect of gap filling alone
    both = per["eligible_peak_picked"] & per["eligible_with_gapfill"]
    r1, r2 = per.loc[both, "rsd_peak_picked_pct"], per.loc[both, "rsd_with_gapfill_pct"]
    rows.append({"mode": "paired_same_features", "n_features_total": len(per), "n_pooled_qc": len(qcs),
                 "min_qc_values_required": need, "n_eligible": int(both.sum()),
                 "median_rsd_pct": f"{r1.median():.2f} vs {r2.median():.2f}",
                 f"n_rsd_le_{limit:g}": f"{int((r1 <= limit).sum())} vs {int((r2 <= limit).sum())}",
                 f"pct_of_eligible_le_{limit:g}": ""})
    return per, pd.DataFrame(rows)


def cal_level(file_id: str) -> str | None:
    import re
    m = re.search(r"ESTD-(A\d)", file_id)
    return m.group(1) if m else None


def is_vs_calibration(ipi: pd.DataFrame, recon_ul: float, spike_ng: float) -> pd.DataFrame:
    rows = []
    for cmp_, g in ipi.groupby("compound"):
        cal = g[(g["role"] == "standard_cal") & g["detected"]].copy()
        cal["level"] = cal["file_id"].map(cal_level)
        cal = cal[cal["level"].notna() & (cal["area"] > 0)]
        cal["ng_per_ul"] = cal["level"].map(ESTD1_NG_PER_UL)
        soil = g[(g["role"] == "sample") & g["detected"] & (g["area"] > 0)]
        row = {"compound": cmp_, "n_cal_points": len(cal), "n_cal_levels": cal["level"].nunique(),
               "n_samples_detected": len(soil), "n_samples": int((g["role"] == "sample").sum())}
        if len(cal["level"].unique()) >= 2 and len(soil):
            x, y = np.log10(cal["ng_per_ul"]), np.log10(cal["area"])
            slope, icpt = np.polyfit(x, y, 1)
            r2 = np.corrcoef(x, y)[0, 1] ** 2
            conc = 10 ** ((np.log10(soil["area"]) - icpt) / slope)
            lo, hi = cal["ng_per_ul"].min(), cal["ng_per_ul"].max()
            found_ng = conc * recon_ul
            row.update({
                "cal_fit": "log10(area) ~ log10(ng/uL)", "slope": round(slope, 3), "r2": round(r2, 4),
                "cal_min_ng_per_ul": lo, "cal_max_ng_per_ul": hi,
                "cal_min_area": cal["area"].min(), "soil_median_area": soil["area"].median(),
                "soil_median_ng_per_ul": conc.median(),
                "n_below_range": int((conc < lo).sum()), "n_above_range": int((conc > hi).sum()),
                "status": "below" if (conc < lo).all() else "above" if (conc > hi).all()
                else "within" if ((conc >= lo) & (conc <= hi)).all() else "partly_below",
                "m_found_median_ng": found_ng.median(),
                "R_app_median_pct": round(float(found_ng.median() / spike_ng * 100), 1),
                "R_app_min_pct": round(float(found_ng.min() / spike_ng * 100), 1),
                "R_app_max_pct": round(float(found_ng.max() / spike_ng * 100), 1),
            })
        else:
            row["status"] = "not_detected_in_samples" if not len(soil) else "too_few_cal_levels"
        rows.append(row)
    return pd.DataFrame(rows)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--results", required=True, help="pipeline results/ directory")
    ap.add_argument("--prefix", required=True, help="e.g. climgrass_pos")
    ap.add_argument("--run-provenance", required=True, help="run_provenance.txt of the pipeline run")
    ap.add_argument("--out", required=True)
    ap.add_argument("--recon-ul", type=float, default=100.0,
                    help="ASSUMED reconstitution volume of soil extracts (uL); confirm")
    ap.add_argument("--spike-ng", type=float, default=10.0)
    a = ap.parse_args(argv)
    res, out = Path(a.results), Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    full_p = res / "mzmine" / f"{a.prefix}.features_full.csv"
    summ_p = res / "qc" / f"{a.prefix}.qc_summary.json"
    inj_p = res / "qc" / f"{a.prefix}.qc_injections.tsv"
    ipi_p = res / "qc" / f"{a.prefix}.qc_is_per_injection.tsv"
    feat_p = res / "qc" / f"{a.prefix}.qc_features.tsv"
    th = json.loads(summ_p.read_text())["thresholds"]

    inj = pd.read_csv(inj_p, sep="\t")
    inj.sort_values("injection_order").to_csv(out / "injections.tsv", sep="\t", index=False)

    full = pd.read_csv(full_p, low_memory=False)
    area = pq.per_file_matrix(full, "area")
    det = pq.detected_matrix(full, area)
    qcs = inj.loc[inj["role"] == "pooled_qc", "file_id"].tolist()
    per, summ = rsd_tables(area, det, qcs, th["qc_min_detect_fraction"], th["qc_feature_rsd_pct"])
    per.index.name = "feature_id"
    per[per["eligible_peak_picked"] | per["eligible_with_gapfill"]].to_csv(out / "rsd_per_feature.tsv", sep="\t")
    summ.to_csv(out / "rsd_summary.tsv", sep="\t", index=False)

    ipi = pd.read_csv(ipi_p, sep="\t")
    med_rt = ipi[ipi["detected"]].groupby("compound")["rt"].median()
    ipi["rt_deviation_min"] = ipi["rt"] - ipi["compound"].map(med_rt)
    ipi = ipi.merge(inj[["file_id", "injection_order", "start_time"]], on="file_id", how="left")
    ipi.sort_values(["compound", "injection_order"]).to_csv(out / "is_per_injection.tsv", sep="\t", index=False)
    blk = ipi[ipi["detected"]].groupby(["compound", "block"]).agg(
        n=("rt", "size"), median_rt=("rt", "median"), median_ppm=("ppm_error", "median"),
        max_abs_ppm=("ppm_error", lambda s: s.abs().max())).reset_index()
    blk.to_csv(out / "is_by_block.tsv", sep="\t", index=False)
    isc = is_vs_calibration(ipi, a.recon_ul, a.spike_ng)
    isc.to_csv(out / "is_vs_calibration.tsv", sep="\t", index=False)

    feat = pd.read_csv(feat_p, sep="\t")
    sd = feat[feat["sample_detected"] > 0]
    br = th["blank_ratio"]
    pd.DataFrame([{
        "n_features_detected_in_samples": len(sd),
        "n_without_blank_signal": int(sd["blank_ratio"].isna().sum()),
        f"n_blank_ratio_gt_{br:g}": int((sd["blank_ratio"] > br).sum()),
        f"pct_blank_ratio_gt_{br:g}": round(float((sd["blank_ratio"] > br).mean() * 100), 1),
        "n_blank_ratio_gt_1": int((sd["blank_ratio"] > 1).sum()),
    }]).to_csv(out / "blank_summary.tsv", sep="\t", index=False)

    shutil.copy(res / "qc" / f"{a.prefix}.qc_drift.tsv", out / "drift.tsv")
    shutil.copy(res / "qc" / f"{a.prefix}.qc_checks.tsv", out / "checks.tsv")
    provenance.write(out, {"features_full": full_p, "qc_summary": summ_p, "qc_injections": inj_p,
                           "qc_is_per_injection": ipi_p, "qc_features": feat_p},
                     {"recon_ul_ASSUMED": a.recon_ul, "spike_ng": a.spike_ng,
                      "estd1_ng_per_ul": ESTD1_NG_PER_UL, "thresholds_from_run": th},
                     notes="Calibration (OE12-1) acquired 27-28 Mar 2023, soils (OE12-3) 30-31 Mar 2023.",
                     pipeline_run=provenance.read_run_provenance(a.run_provenance))


if __name__ == "__main__":
    main()
