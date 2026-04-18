#!/usr/bin/env python3
"""Build the SIMPER fingerprint atlas from the consensus-aligned feature table.

STATUS: v0.1 stub. Real implementation deferred to v0.2.

For v0.1 demo runs the SIMPER fingerprint atlas is consumed pre-built from
Zenodo (see ``bin/fetch_atlas.py``) - it is the canonical required artefact
downstream decomposition reads. Invoking this wrapper still exits non-zero
with a clear stderr message so Nextflow surfaces a process failure if a
caller accidentally flips ``mode = 'train'`` in v0.1.

To implement in v0.2, port the SIMPER top-N per-phylum selection code from
``soilmass-viewer/analysis/analysis-15/05_simper/`` keeping the same
direction (UP/DOWN) and fold-change outputs downstream scripts expect.

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
    p.add_argument("--top-n", type=int, default=50)
    p.add_argument("--min-batches", type=int, default=3)
    p.add_argument(
        "--output-simper-fp-atlas",
        required=True,
        help="Destination parquet for the SIMPER fingerprint atlas.",
    )
    p.parse_args()

    print(
        "build_simper_fp_atlas.py is a v0.1 stub. SIMPER fingerprint training "
        "is deferred to v0.2. For v0.1 demo runs, fetch the pre-built atlas "
        "via fetch_atlas.py.",
        file=sys.stderr,
    )
    return NOT_IMPLEMENTED_EXIT


if __name__ == "__main__":
    sys.exit(main() or 0)
