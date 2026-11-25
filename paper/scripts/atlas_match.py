#!/usr/bin/env python3
"""Issue #6: minimal stage 7 - match ClimGrass MS/MS spectra to the Atlas reference library.

Criteria (manuscript 2.8): precursor within 10 ppm, modified cosine >= 0.7 with
fragment tolerance 0.02 Da, at least 4 matched fragments (matchms ModifiedCosineGreedy,
top 50 peaks per spectrum, as bin/simper_match.py).

Library: reference MS/MS of the release reprocessing (one spectrum per reference
feature), limited to features detected in at least one reference sample, so soil
injections in the reprocessing cannot match themselves.

Writes:
  matches.tsv            every passing (soil feature, reference feature) pair
  best_match.tsv         best match per soil feature, with its phyla
  coverage.tsv           soil features with MS/MS, with a match, share
  phylum_detection.tsv   per matched reference feature: samples detected / samples per phylum
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import provenance  # noqa: E402


def candidates(query_mz: float, lib_mz: np.ndarray, ppm: float) -> np.ndarray:
    """Indices of library precursors within ppm of query_mz; lib_mz must be sorted."""
    tol = query_mz * ppm * 1e-6
    lo, hi = np.searchsorted(lib_mz, query_mz - tol, "left"), np.searchsorted(lib_mz, query_mz + tol, "right")
    return np.arange(lo, hi)


def detection_by_phylum(presence: pd.DataFrame, features: list) -> pd.DataFrame:
    """presence: feature_id, unit, phylum. Rate = units with the feature / units of that phylum."""
    per_phylum = presence.drop_duplicates(["unit", "phylum"]).groupby("phylum")["unit"].nunique()
    sub = presence[presence["feature_id"].isin(features)]
    det = sub.drop_duplicates(["feature_id", "unit"]).groupby(["feature_id", "phylum"])["unit"].nunique()
    out = det.rename("n_detected").reset_index()
    out["n_units_phylum"] = out["phylum"].map(per_phylum)
    out["detection_rate"] = out["n_detected"] / out["n_units_phylum"]
    return out


def _load(path, keep=None):
    from matchms.filtering import default_filters, normalize_intensities, select_by_relative_intensity, reduce_to_number_of_peaks
    from matchms.importing import load_from_mgf
    out = []
    for s in load_from_mgf(str(path)):
        fid = s.get("feature_id") or s.get("scans")
        if keep is not None and str(fid) not in keep:
            continue
        s = default_filters(s)
        s = normalize_intensities(s)
        s = reduce_to_number_of_peaks(s, n_max=50)
        if s is None or s.get("precursor_mz") is None or len(s.peaks.mz) == 0:
            continue
        s.set("feature_id", str(fid))
        out.append(s)
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--soil-mgf", required=True)
    ap.add_argument("--library-mgf", required=True)
    ap.add_argument("--presence", required=True, help="feature_id,unit,phylum (release a_presence_long.csv)")
    ap.add_argument("--run-provenance", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--ppm", type=float, default=10.0)
    ap.add_argument("--tol", type=float, default=0.02)
    ap.add_argument("--min-cos", type=float, default=0.7)
    ap.add_argument("--min-peaks", type=int, default=4)
    a = ap.parse_args(argv)
    from matchms.similarity import ModifiedCosineGreedy
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)

    presence = pd.read_csv(a.presence, dtype={"feature_id": str})
    ref_ids = set(presence["feature_id"])
    soil = _load(a.soil_mgf)
    lib = _load(a.library_mgf, keep=ref_ids)
    lib.sort(key=lambda s: s.get("precursor_mz"))
    lib_mz = np.array([s.get("precursor_mz") for s in lib])
    mc = ModifiedCosineGreedy(tolerance=a.tol)

    rows = []
    for q in soil:
        for i in candidates(q.get("precursor_mz"), lib_mz, a.ppm):
            r = mc.pair(lib[i], q)
            score, n = float(r["score"]), int(r["matches"])
            if score >= a.min_cos and n >= a.min_peaks:
                rows.append({"soil_feature_id": q.get("feature_id"), "soil_mz": q.get("precursor_mz"),
                             "ref_feature_id": lib[i].get("feature_id"), "ref_mz": lib_mz[i],
                             "ppm": (lib_mz[i] - q.get("precursor_mz")) / q.get("precursor_mz") * 1e6,
                             "cosine": round(score, 4), "matched_peaks": n})
    m = pd.DataFrame(rows, columns=["soil_feature_id", "soil_mz", "ref_feature_id", "ref_mz", "ppm", "cosine",
                                    "matched_peaks"])
    m.to_csv(out / "matches.tsv", sep="\t", index=False)
    best = m.sort_values(["cosine", "matched_peaks"], ascending=False).drop_duplicates("soil_feature_id")
    phy = presence.drop_duplicates(["feature_id", "phylum"]).groupby("feature_id")["phylum"].agg(
        lambda s: ";".join(sorted(s)))
    best = best.assign(ref_phyla=best["ref_feature_id"].map(phy),
                       n_ref_phyla=best["ref_feature_id"].map(phy).str.count(";") + 1)
    best.to_csv(out / "best_match.tsv", sep="\t", index=False)
    detection_by_phylum(presence, best["ref_feature_id"].unique().tolist()).to_csv(
        out / "phylum_detection.tsv", sep="\t", index=False)
    pd.DataFrame([{"soil_features_with_ms2": len(soil), "library_spectra": len(lib),
                   "soil_features_matched": best["soil_feature_id"].nunique(),
                   "pct_matched": round(best["soil_feature_id"].nunique() / max(len(soil), 1) * 100, 1),
                   "reference_features_hit": best["ref_feature_id"].nunique(),
                   "n_phyla_in_presence": presence["phylum"].nunique()}]).to_csv(out / "coverage.tsv", sep="\t",
                                                                                     index=False)
    provenance.write(out, {"soil_mgf": a.soil_mgf, "library_mgf": a.library_mgf, "presence": a.presence},
                     {"ppm": a.ppm, "fragment_tol_da": a.tol, "min_cosine": a.min_cos, "min_matched_peaks": a.min_peaks,
                      "top_n_peaks": 50, "similarity": "matchms ModifiedCosineGreedy"},
                     notes="Library = release reprocessing sla_reprocess/t3 (mzmine_full/full_ann20_fbmn.mgf) "
                           "limited to features in analysis_pos/a_presence_long.csv.",
                     pipeline_run=provenance.read_run_provenance(a.run_provenance))


if __name__ == "__main__":
    main()
