#!/usr/bin/env python3
"""Top-N contributing features per phylum (v0.1 diagnostic).

Scope note: this is the v0.1 minimal version. It produces a ranked list of
the dominant features per phylum from the SIMPER fingerprint atlas, paired
with each feature's mean intensity across soil samples from the corrected
intensity matrix. Intended for eyeball inspection before every release —
if a single low-mass artefact is topping a phylum's list, that's the signal
that the v0.2 "RIE amplification detector" (described in analysis-19 notes)
needs to run.

Inputs:
  --corrected-intensity  Parquet, features x samples (feature_id as index).
                         Typically the output of apply_rie_correction.
  --simper-atlas         SIMPER fingerprint atlas parquet. Must have columns
                         feature_id, phylum, simper_rank (ascending = top
                         contributor). The optional `feature_class` column is
                         emitted if present.
  --top-n                How many features per phylum to report (default 10).
  --output               Output TSV with columns:
                         phylum, rank, feature_id,
                         mean_intensity_across_samples,
                         feature_class_if_annotated
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--corrected-intensity", required=True, type=Path,
                   help="Corrected intensity parquet (features x samples)")
    p.add_argument("--simper-atlas", required=True, type=Path,
                   help="SIMPER fingerprint atlas parquet "
                        "(columns: feature_id, phylum, simper_rank)")
    p.add_argument("--top-n", type=int, default=10,
                   help="Number of top features per phylum to report (default: 10)")
    p.add_argument("--output", required=True, type=Path,
                   help="Output TSV")
    args = p.parse_args(argv)

    intensity = pd.read_parquet(args.corrected_intensity)
    simper = pd.read_parquet(args.simper_atlas)

    required_simper = {"feature_id", "phylum", "simper_rank"}
    if not required_simper.issubset(simper.columns):
        sys.exit(
            f"simper atlas missing columns; got {list(simper.columns)}, "
            f"need at least {sorted(required_simper)}"
        )

    # Mean intensity per feature across samples. Features missing from the
    # intensity matrix get NaN — callers see that as a signal that the atlas
    # carries features the current corrected matrix can't speak to.
    mean_intensity = intensity.mean(axis=1)

    has_feature_class = "feature_class" in simper.columns

    rows: list[dict] = []
    for phylum, sub in simper.groupby("phylum"):
        # Ascending simper_rank = most-enriched first. Use stable sort so
        # ties land in whatever deterministic order the atlas provided.
        top = sub.sort_values("simper_rank", kind="mergesort").head(args.top_n)
        for _, r in top.iterrows():
            fid = r["feature_id"]
            row = {
                "phylum": phylum,
                "rank": int(r["simper_rank"]),
                "feature_id": fid,
                "mean_intensity_across_samples": (
                    float(mean_intensity[fid]) if fid in mean_intensity.index else float("nan")
                ),
            }
            if has_feature_class:
                row["feature_class_if_annotated"] = r["feature_class"]
            else:
                row["feature_class_if_annotated"] = ""
            rows.append(row)

    out_df = pd.DataFrame(rows, columns=[
        "phylum", "rank", "feature_id",
        "mean_intensity_across_samples",
        "feature_class_if_annotated",
    ])
    args.output.parent.mkdir(parents=True, exist_ok=True)
    out_df.to_csv(args.output, sep="\t", index=False)
    return 0


if __name__ == "__main__":
    sys.exit(main())
