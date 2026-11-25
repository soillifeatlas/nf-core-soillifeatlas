#!/usr/bin/env python3
"""Issues #6 and #8: mirror-plot spectra and per-phylum detection for chosen soil-to-Atlas
matches, and Atlas occurrence of the drought-responsive lipids.

  mirror_<soil_feature_id>.tsv   soil and reference peaks (m/z, relative intensity, side)
  examples.tsv                   the chosen matches with cosine, fragments, phyla
  examples_phylum_detection.tsv  detection rate per phylum for the chosen reference features
  responsive_atlas_occurrence.tsv  responsive lipids (08_application) with their best Atlas match
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import provenance  # noqa: E402


def read_mgf(path: str, wanted: set[str]) -> dict[str, pd.DataFrame]:
    """Peaks of the spectra whose FEATURE_ID is in wanted (streams the file)."""
    out, cur, peaks, keep = {}, None, [], False
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line == "BEGIN IONS":
                cur, peaks, keep = None, [], False
            elif line.startswith("FEATURE_ID="):
                cur = line.split("=", 1)[1]
                keep = cur in wanted
            elif line == "END IONS":
                if keep and cur is not None:
                    d = pd.DataFrame(peaks, columns=["mz", "intensity"])
                    d["rel_intensity"] = d["intensity"] / d["intensity"].max() * 100
                    out[cur] = d
            elif keep and line and line[0].isdigit():
                mz, it = line.split()[:2]
                peaks.append((float(mz), float(it)))
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--match-dir", required=True, help="paper/06_atlas_match/pos")
    ap.add_argument("--soil-mgf", required=True)
    ap.add_argument("--library-mgf", required=True)
    ap.add_argument("--annotations", required=True, help="paper/02_annotation/pos/annotations.tsv")
    ap.add_argument("--responsive", required=True, help="paper/08_application/pos/responsive.tsv")
    ap.add_argument("--examples", required=True, help="comma-separated soil feature ids, broad first")
    ap.add_argument("--out", required=True)
    a = ap.parse_args(argv)
    md, out = Path(a.match_dir), Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    best = pd.read_csv(md / "best_match.tsv", sep="\t", dtype={"soil_feature_id": str, "ref_feature_id": str})
    det = pd.read_csv(md / "phylum_detection.tsv", sep="\t", dtype={"feature_id": str})
    ann = pd.read_csv(a.annotations, sep="\t", dtype={"feature_id": str})
    ids = a.examples.split(",")
    ex = best[best["soil_feature_id"].isin(ids)].merge(
        ann[["feature_id", "lipid_name", "confidence", "rt"]], left_on="soil_feature_id", right_on="feature_id")
    ex["role"] = ex["soil_feature_id"].map({ids[0]: "broad", **{i: "restricted" for i in ids[1:]}})
    ex.drop(columns="feature_id").to_csv(out / "examples.tsv", sep="\t", index=False)
    det[det["feature_id"].isin(ex["ref_feature_id"])].to_csv(out / "examples_phylum_detection.tsv", sep="\t",
                                                             index=False)
    soil = read_mgf(a.soil_mgf, set(ids))
    lib = read_mgf(a.library_mgf, set(ex["ref_feature_id"]))
    for _, r in ex.iterrows():
        s = soil[r["soil_feature_id"]].assign(side="soil")
        ref = lib[r["ref_feature_id"]].assign(side="atlas")
        pd.concat([s, ref]).to_csv(out / f"mirror_{r['soil_feature_id']}.tsv", sep="\t", index=False)

    resp = pd.read_csv(a.responsive, sep="\t", dtype={"feature_id": str})
    occ = resp.merge(best, left_on="feature_id", right_on="soil_feature_id", how="left")
    occ[["feature_id", "lipid_name", "lipid_class", "confidence", "quant_route", "log2fc_drought_vs_ambient",
         "q_drought", "ref_feature_id", "cosine", "matched_peaks", "n_ref_phyla", "ref_phyla"]].to_csv(
        out / "responsive_atlas_occurrence.tsv", sep="\t", index=False)
    provenance.write(out, {"best_match": md / "best_match.tsv", "phylum_detection": md / "phylum_detection.tsv",
                           "soil_mgf": a.soil_mgf, "library_mgf": a.library_mgf, "annotations": a.annotations,
                           "responsive": a.responsive},
                     {"examples": ids, "selection": "Gold, unflagged; broad = most phyla, restricted = one phylum, "
                                                    "highest cosine (proposed, to confirm)"})


if __name__ == "__main__":
    main()
