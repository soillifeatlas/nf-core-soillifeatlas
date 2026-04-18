#!/usr/bin/env python3
"""Cross-batch alignment via 3-stage LOESS pipeline (anchor -> LOESS -> mapping).

STATUS: v0.1 stub. Real implementation deferred to v0.2.

For v0.1 demo runs the atlas is consumed pre-aligned from Zenodo (see
``bin/fetch_atlas.py``), so the train workflow that would call this script
is scaffolded but not demo-run. Invoking this wrapper still exits non-zero
with a clear stderr message so Nextflow surfaces a process failure if a
caller accidentally flips ``mode = 'train'`` in v0.1.

To implement in v0.2, port the 3-stage alignment code from
``soilmass-viewer/analysis/analysis-15/03_alignment/`` using these params:

  * anchor:  5 ppm m/z, ±2.0 min RT, >=3/6 batches, RT drift < 1.0 min
  * LOESS:   frac 0.3, locally linear, degree 1
  * mapping: ±0.5 min RT, 5 ppm, RT window 1.5 - 25.0 min

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
        "--inputs",
        nargs="+",
        required=True,
        help="Per-batch feature-table CSVs to be aligned.",
    )
    p.add_argument("--anchor-ppm", type=float, default=5.0)
    p.add_argument("--anchor-rt-window", type=float, default=2.0)
    p.add_argument("--min-batches", type=int, default=3)
    p.add_argument("--loess-frac", type=float, default=0.3)
    p.add_argument("--map-rt-window", type=float, default=0.5)
    p.add_argument("--map-ppm", type=float, default=5.0)
    p.add_argument("--rt-min", type=float, default=1.5)
    p.add_argument("--rt-max", type=float, default=25.0)
    p.add_argument(
        "--output",
        required=True,
        help="Destination parquet for the consensus-aligned feature table.",
    )
    p.parse_args()

    print(
        "align_loess.py is a v0.1 stub. Cross-batch alignment is deferred to v0.2. "
        "For v0.1 demo runs, atlas data is pre-aligned - fetch it via fetch_atlas.py.",
        file=sys.stderr,
    )
    return NOT_IMPLEMENTED_EXIT


if __name__ == "__main__":
    sys.exit(main() or 0)
