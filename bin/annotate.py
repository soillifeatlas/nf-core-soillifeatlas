#!/usr/bin/env python3
"""Stage 4: merge lipid annotations into one table with confidence levels.

Inputs: the MZmine full feature table of one polarity (always) and, as
plug-ins, any of: a LipidSearch export, an MS-DIAL alignment export, spectral
library matches (GNPS library results format) and a SIRIUS/CANOPUS summary.
Outputs, all prefixed with --prefix:

  <prefix>.annotation_sources.tsv   every candidate from every source
  <prefix>.annotations.tsv          one row per annotated feature: name,
                                    confidence (Gold/Silver/Bronze), sources,
                                    quantification route (Q1-Q3)
  <prefix>.annotation_summary.json  counts per confidence, source and class
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import pandas as pd

from soillifeatlas.annotate import (
    check_rt,
    from_library,
    from_lipidsearch,
    from_msdial,
    from_mzmine,
    from_sirius,
    merge,
    read_table,
    summarise,
)


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--features", required=True, help="MZmine full feature table (CSV)")
    p.add_argument("--polarity", required=True, choices=["pos", "neg"])
    p.add_argument("--prefix", required=True)
    p.add_argument("--lipidsearch", default="", help="LipidSearch export (CSV, TSV or XLSX)")
    p.add_argument("--msdial", default="", help="MS-DIAL alignment export (TSV)")
    p.add_argument("--library", default="", help="Spectral library matches, GNPS results format")
    p.add_argument("--library-same-run", action="store_true",
                   help="Library #Scan# values are this run's feature ids")
    p.add_argument("--sirius", default="", help="SIRIUS/CANOPUS summary (TSV, CSV or XLSX)")
    p.add_argument("--is-targets", default="", help="Internal-standard table; its classes allow Q1")
    p.add_argument("--mz-ppm", type=float, default=5.0)
    p.add_argument("--rt-window", type=float, default=0.2)
    p.add_argument("--library-min-score", type=float, default=0.7)
    p.add_argument("--library-min-peaks", type=int, default=4)
    return p.parse_args(argv)


def main(argv=None) -> int:
    a = parse_args(argv)
    # a plug-in that is a placeholder or an empty file counts as not supplied
    for k in ("lipidsearch", "msdial", "library", "sirius"):
        p = getattr(a, k)
        if p and (Path(p).name.startswith("NO_FILE") or Path(p).stat().st_size == 0):
            setattr(a, k, "")
    keep = ["id", "mz", "rt", "fragment_scans", "lipid_annotations:lipid_annotations",
            "lipid_annotations:ion_adduct", "lipid_annotations:explained_intensity_percent"]
    head = pd.read_csv(a.features, nrows=0).columns
    full = pd.read_csv(a.features, usecols=[c for c in keep if c in head], low_memory=False)
    features = full[["id", "mz", "rt"] + (["fragment_scans"] if "fragment_scans" in full else [])]

    parts, inputs = [from_mzmine(full)], {"MZmine lipid rules": int(full.get(keep[4], pd.Series()).notna().sum())}
    if a.lipidsearch:
        t = read_table(a.lipidsearch, a.polarity)
        parts.append(from_lipidsearch(t, features, a.mz_ppm, a.rt_window))
        inputs["LipidSearch"] = len(t)
    if a.msdial:
        t = read_table(a.msdial, a.polarity)
        parts.append(from_msdial(t, features, a.mz_ppm, a.rt_window))
        inputs["MS-DIAL"] = len(t)
    if a.library:
        t = read_table(a.library, a.polarity)
        parts.append(from_library(t, features, a.mz_ppm, a.rt_window, a.library_min_score,
                                  a.library_min_peaks, a.library_same_run))
        inputs["spectral library"] = len(t)
    if a.sirius:
        t = read_table(a.sirius, a.polarity)
        parts.append(from_sirius(t, features, a.mz_ppm, a.rt_window))
        inputs["SIRIUS/CANOPUS"] = len(t)
    cands = pd.concat([p for p in parts if len(p)], ignore_index=True) if any(len(p) for p in parts) \
        else parts[0]

    is_classes = set()
    if a.is_targets:
        tg = pd.read_csv(a.is_targets)
        if "lipid_class" in tg.columns:
            is_classes = set(tg["lipid_class"].dropna())

    merged = check_rt(merge(cands, features, is_classes))
    summary = summarise(merged, cands, features, inputs)
    summary.update({"prefix": a.prefix, "polarity": a.polarity, "mz_ppm": a.mz_ppm, "rt_window": a.rt_window})

    out = Path(".")
    cands.to_csv(out / f"{a.prefix}.annotation_sources.tsv", sep="\t", index=False)
    merged.to_csv(out / f"{a.prefix}.annotations.tsv", sep="\t", index=False)
    (out / f"{a.prefix}.annotation_summary.json").write_text(json.dumps(summary, indent=2, default=str))
    c = summary["confidence"]
    print(f"Annotation {a.prefix}: {summary['annotated_features']} features annotated "
          f"(Gold {c['Gold']}, Silver {c['Silver']}, Bronze {c['Bronze']})")
    return 0


if __name__ == "__main__":
    sys.exit(main())
