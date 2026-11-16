#!/usr/bin/env python3
"""Build the per-phylum biomarker atlas from the consensus-aligned feature table.

STATUS: v0.1 stub. Real implementation deferred to v0.2.

For v0.1 demo runs the biomarker atlas is consumed pre-built from Zenodo
(see ``bin/fetch_atlas.py``). Invoking this wrapper still exits non-zero
with a clear stderr message so Nextflow surfaces a process failure if a
caller accidentally flips ``mode = 'train'`` in v0.1.

To implement in v0.2, port the Indval+composite biomarker code from
``soilmass-viewer/analysis/analysis-15/04_biomarkers/`` honouring the
design's thresholds (Indval IndA >= threshold across >= min-batches).

The CLI surface below mirrors the expected v0.2 contract so Nextflow
modules can be written against it today.
"""
from __future__ import annotations

import argparse
import sys


NOT_IMPLEMENTED_EXIT = 2


def main() -> int:
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument(
        "--consensus-aligned-table",
        required=True,
        help="Input parquet from align_loess.py",
    )
    p.add_argument(
        "--sample-metadata",
        required=True,
        help="Per-sample metadata (batch, kingdom, phylum, etc.).",
    )
    p.add_argument(
        "--output-biomarker-atlas",
        required=True,
        help="Destination parquet for the per-phylum biomarker atlas.",
    )
    p.add_argument(
        "--composite-weights-json",
        default=None,
        help="Optional composite scoring weights.",
    )
    p.add_argument("--indval-threshold", type=float, default=0.7)
    p.add_argument("--indval-min-batches", type=int, default=3)
    p.parse_args()

    print(
        "build_biomarker_atlas.py is a v0.1 stub. Biomarker-atlas training is "
        "deferred to v0.2. For v0.1 demo runs, fetch the pre-built atlas via "
        "fetch_atlas.py.",
        file=sys.stderr,
    )
    return NOT_IMPLEMENTED_EXIT


if __name__ == "__main__":
    sys.exit(main() or 0)
