#!/usr/bin/env python3
"""Build assets/equisplash_targets.csv from assets/equisplash13_targets_source.json.

The source JSON lists the 13 EquiSPLASH components with formula, one or two
"primary" ions per polarity, their m/z and the spiked amount (ng). This script
derives each component's neutral monoisotopic mass from its first listed ion
and expands it to the common LC-MS adducts, so the QC step can still find a
standard that ionised as a different adduct (for example [M+Na]+ instead of
[M+H]+, which happens for PE, PS and LPE in the ClimGrass data).

Ions listed in the source are marked priority=primary; derived ones secondary.

Run from the repository root:
    python3 scripts/make_equisplash_targets.py
"""
from __future__ import annotations

import csv
import json
from pathlib import Path

# Adduct mass shifts (Da) relative to the neutral monoisotopic mass M,
# electron mass included.
ADDUCTS = {
    "pos": {
        "[M+H]+": 1.007276,
        "[M+NH4]+": 18.033823,
        "[M+Na]+": 22.989218,
        "[M+H-H2O]+": -17.003289,
    },
    "neg": {
        "[M-H]-": -1.007276,
        "[M+HCOO]-": 44.998201,
    },
}

ROOT = Path(__file__).resolve().parent.parent
SRC = ROOT / "assets" / "equisplash13_targets_source.json"
OUT = ROOT / "assets" / "equisplash_targets.csv"


def short_name(name: str) -> str:
    """'IS PC 15:0-18:1(d7)' -> 'PC_15_0_18_1_d7' (matches the pipeline's IS ids)."""
    core = name.removeprefix("IS ").strip()
    for ch in " :-()":
        core = core.replace(ch, "_")
    while "__" in core:
        core = core.replace("__", "_")
    return core.strip("_")


def main() -> None:
    src = json.loads(SRC.read_text())
    compounds: dict[str, dict] = {}
    primary: set[tuple[str, str, str]] = set()
    for entry in src:
        cid = short_name(entry["name"])
        pol = entry["polarity"].lower()
        shift = ADDUCTS[pol][entry["ion"]]
        neutral = entry["mz"] - shift
        compounds.setdefault(
            cid,
            {
                "compound": cid,
                "name": entry["name"].removeprefix("IS ").strip(),
                "lipid_class": entry["name"].removeprefix("IS ").split()[0],
                "formula": entry["formula"],
                "neutral_mass": neutral,
                "spike_ng": entry.get("spike_ng", ""),
            },
        )
        primary.add((cid, pol, entry["ion"]))

    rows = []
    for cid, c in sorted(compounds.items()):
        for pol, adducts in ADDUCTS.items():
            for ion, shift in adducts.items():
                rows.append(
                    {
                        "compound": cid,
                        "name": c["name"],
                        "lipid_class": c["lipid_class"],
                        "formula": c["formula"],
                        "polarity": pol,
                        "adduct": ion,
                        "mz": f"{c['neutral_mass'] + shift:.5f}",
                        "priority": "primary" if (cid, pol, ion) in primary else "secondary",
                        "spike_ng": c["spike_ng"],
                    }
                )
    with OUT.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    print(f"wrote {len(rows)} targets for {len(compounds)} compounds to {OUT}")


if __name__ == "__main__":
    main()
