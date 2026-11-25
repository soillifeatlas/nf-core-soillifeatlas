#!/usr/bin/env python3
"""Write copies of mzML files that keep every MS1 scan and only one MS2 activation type.

    keep_ms2_activation.py --activation hcd --out-dir filtered/ data/*.mzML

The copies keep their file names. A table of scan counts per file is written to
<out-dir>/activation_filter.tsv, and a file that keeps no fragment scan is reported.
See soillifeatlas.mzml_filter for why.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

from soillifeatlas.mzml_filter import keep_activation


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--activation', required=True, choices=['hcd', 'cid'])
    ap.add_argument('--out-dir', required=True, type=Path)
    ap.add_argument('files', nargs='+', type=Path)
    a = ap.parse_args(argv)
    a.out_dir.mkdir(parents=True, exist_ok=True)
    rows = ['file\tms1\tkept_msn\tdropped_msn']
    empty = []
    for f in a.files:
        out, n = keep_activation(f.read_text(encoding='utf-8'), a.activation)
        (a.out_dir / f.name).write_text(out, encoding='utf-8')
        rows.append(f"{f.name}\t{n['ms1']}\t{n['kept_msn']}\t{n['dropped_msn']}")
        if n['kept_msn'] == 0:
            empty.append(f.name)
    (a.out_dir / 'activation_filter.tsv').write_text('\n'.join(rows) + '\n', encoding='utf-8')
    print(f"{len(a.files)} files written to {a.out_dir}")
    if empty:
        print(f"no {a.activation} fragment scans in: {', '.join(empty)}", file=sys.stderr)
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
