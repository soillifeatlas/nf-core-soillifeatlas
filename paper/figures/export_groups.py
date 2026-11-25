#!/usr/bin/env python
"""Export what the figures need from an Atlas release.json into paper/.

    python paper/figures/export_groups.py <path/to/release.json>

Writes, under paper/figures/inputs/:
  atlas_groups.csv          phylum -> organism group, reference samples per phylum
  atlas_release_counts.csv  headline counts of the release (key, value)
  atlas_feature_layers.csv  feature layers from consensus features to annotated Atlas lipids
  atlas_groups_provenance.yaml  release label, source hash, export time
The figure scripts read only these files (inside paper/), never the release itself.
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
    with open(OUT / "atlas_release_counts.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["key", "value"])
        w.writerows(r["counts"].items())
        kf = r.get("key_finding", {})
        for k in ("exclusive_features", "features_considered"):
            if k in kf:
                w.writerow([f"key_finding_{k}", kf[k]])
    with open(OUT / "atlas_feature_layers.csv", "w", newline="") as fh:
        cols = ["id", "label", "pos", "neg", "derived_from", "used_for"]
        w = csv.DictWriter(fh, fieldnames=cols, extrasaction="ignore")
        w.writeheader()
        w.writerows({c: layer.get(c) for c in cols} for layer in r["feature_layers"])
    (OUT / "atlas_groups_provenance.yaml").write_text(yaml.safe_dump({
        "exported": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),
        "command": "python paper/figures/export_groups.py <release.json>",
        "source_file": "SoilLifeAtlas/release/release.json",
        "source_sha256": hashlib.sha256(raw).hexdigest(),
        "release": r["release"], "taxonomy": r["taxonomy"],
        "n_groups": len(r["groups"]), "n_phyla": len(rows),
        "key_finding": r.get("key_finding"),
        "files": ["atlas_groups.csv", "atlas_release_counts.csv", "atlas_feature_layers.csv"],
    }, sort_keys=False, allow_unicode=True))
    print(f"wrote {OUT / 'atlas_groups.csv'} ({len(rows)} phyla)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1]))
