#!/usr/bin/env python3
"""Mann-Whitney U test per (phylum x contrast) for treatment-effect checks.

Ported (not wrapped) — the framework library does not have an
equivalent; the test is a straightforward two-sample rank-sum so we
implement it here against scipy.

For each (phylum, contrast_factor) pair, partitions samples into
`Yes`/`No` groups based on the metadata column named by `--contrast`,
then runs a two-sided Mann-Whitney U on the two proportion vectors.
Direction is derived from the median difference:
  - "+"  if median(Yes) > median(No)
  - "-"  if median(Yes) < median(No)
  - "0"  if the medians are equal (exact tie)

Inputs:
  --composition-phylum  Parquet in long format (sample_id, phylum,
                        proportion_pct).
  --sample-metadata     TSV with at minimum a `sample_id` column plus one
                        column per contrast factor. Values in the contrast
                        columns must be exactly "Yes"/"No".
  --contrast            Space-separated list of metadata-column names to
                        test (e.g. `--contrast drought warming`).
  --output              Output TSV with columns:
                        phylum, contrast, direction, p_value,
                        median_yes, median_no, n_yes, n_no
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu


def _direction(median_yes: float, median_no: float) -> str:
    """+/-/0 direction encoding for median(Yes) vs median(No) comparison."""
    if median_yes > median_no:
        return "+"
    if median_yes < median_no:
        return "-"
    return "0"


def _mwu_for_contrast(
    comp: pd.DataFrame,
    meta: pd.DataFrame,
    contrast: str,
) -> list[dict]:
    """Run Mann-Whitney U for every phylum under the given contrast column.

    Returns one row per phylum with the stats flattened into a dict.
    """
    if contrast not in meta.columns:
        sys.exit(
            f"contrast factor {contrast!r} not found in sample metadata columns "
            f"(available: {sorted(meta.columns)})"
        )
    yes_samples = set(meta.loc[meta[contrast] == "Yes", "sample_id"])
    no_samples = set(meta.loc[meta[contrast] == "No", "sample_id"])

    rows = []
    for phylum, sub in comp.groupby("phylum"):
        pvt = sub.set_index("sample_id")["proportion_pct"]
        group_yes = pvt.reindex([s for s in pvt.index if s in yes_samples]).dropna().values
        group_no = pvt.reindex([s for s in pvt.index if s in no_samples]).dropna().values

        n_yes = int(len(group_yes))
        n_no = int(len(group_no))

        if n_yes == 0 or n_no == 0:
            # Undefined under this contrast — emit NaNs, skip the test.
            rows.append({
                "phylum": phylum,
                "contrast": contrast,
                "direction": "0",
                "p_value": float("nan"),
                "median_yes": float("nan") if n_yes == 0 else float(np.median(group_yes)),
                "median_no":  float("nan") if n_no == 0 else float(np.median(group_no)),
                "n_yes": n_yes,
                "n_no": n_no,
            })
            continue

        median_yes = float(np.median(group_yes))
        median_no = float(np.median(group_no))

        try:
            _stat, pval = mannwhitneyu(group_yes, group_no, alternative="two-sided")
            pval = float(pval)
        except ValueError:
            # scipy raises if both groups are constant and identical — a
            # legitimate "no separation" signal, not an error condition.
            pval = float("nan")

        rows.append({
            "phylum": phylum,
            "contrast": contrast,
            "direction": _direction(median_yes, median_no),
            "p_value": pval,
            "median_yes": median_yes,
            "median_no": median_no,
            "n_yes": n_yes,
            "n_no": n_no,
        })
    return rows


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--composition-phylum", required=True, type=Path,
                   help="Long-format phylum composition parquet "
                        "(sample_id, phylum, proportion_pct)")
    p.add_argument("--sample-metadata", required=True, type=Path,
                   help="TSV with sample_id + one column per contrast factor "
                        "(values must be Yes/No)")
    p.add_argument("--contrast", required=True, nargs="+",
                   help="Space-separated list of metadata-column names "
                        "to use as contrast factors (e.g. drought warming)")
    p.add_argument("--output", required=True, type=Path,
                   help="Output TSV (one row per phylum x contrast pair)")
    args = p.parse_args(argv)

    comp = pd.read_parquet(args.composition_phylum)
    meta = pd.read_csv(args.sample_metadata, sep="\t")

    required_comp = {"sample_id", "phylum", "proportion_pct"}
    if not required_comp.issubset(comp.columns):
        sys.exit(
            f"composition_phylum missing columns; got {list(comp.columns)}, "
            f"need {sorted(required_comp)}"
        )
    if "sample_id" not in meta.columns:
        sys.exit(
            f"sample metadata must have a 'sample_id' column; got {list(meta.columns)}"
        )

    all_rows: list[dict] = []
    for contrast in args.contrast:
        all_rows.extend(_mwu_for_contrast(comp, meta, contrast))

    out_df = pd.DataFrame(all_rows, columns=[
        "phylum", "contrast", "direction", "p_value",
        "median_yes", "median_no", "n_yes", "n_no",
    ])
    args.output.parent.mkdir(parents=True, exist_ok=True)
    out_df.to_csv(args.output, sep="\t", index=False)
    return 0


if __name__ == "__main__":
    sys.exit(main())
