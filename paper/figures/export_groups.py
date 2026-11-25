#!/usr/bin/env python
"""Export the Atlas phylum -> organism-group mapping from a release.json into paper/.

    python paper/figures/export_groups.py <path/to/release.json>

Writes paper/figures/inputs/atlas_groups.csv and atlas_groups_provenance.yaml. The figure
scripts read only that CSV (inside paper/), never the release itself.
"""
from __future__ import annotations

import csv
import datetime as dt
import hashlib
import json
import sys
from pathlib import Path

import yaml

OUT = Path(__file__).resolve().parent / "inputs"


def main(release_json: str) -> int:
    src = Path(release_json).expanduser().resolve()
    raw = src.read_bytes()
    r = json.loads(raw)
    OUT.mkdir(exist_ok=True)
    rows = []
    for g in r["groups"]:
        for p in g["phyla"]:
            rows.append({"group_id": g["id"], "group": g["name"], "phylum": p["name"],
                         "samples_pos": p.get("samples_pos"), "samples_neg": p.get("samples_neg"),
                         "analysed": p.get("analysed")})
    with open(OUT / "atlas_groups.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    (OUT / "atlas_groups_provenance.yaml").write_text(yaml.safe_dump({
        "exported": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),
        "command": "python paper/figures/export_groups.py <release.json>",
        "source_file": "SoilLifeAtlas/release/release.json",
        "source_sha256": hashlib.sha256(raw).hexdigest(),
        "release": r["release"], "taxonomy": r["taxonomy"],
        "n_groups": len(r["groups"]), "n_phyla": len(rows),
    }, sort_keys=False, allow_unicode=True))
    print(f"wrote {OUT / 'atlas_groups.csv'} ({len(rows)} phyla)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1]))
