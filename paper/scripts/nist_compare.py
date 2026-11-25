#!/usr/bin/env python3
"""Issue #5: compare stage-5 estimates for NIST SRM 1950 with the interlaboratory consensus.

Consensus: Bowden et al. (2017) J Lipid Res 58:2275 / NISTIR 8185, Tables 2-6 (lipids
measured by >= 5 laboratories), final MEDM location in nmol/mL.

  --parse   NISTIR text (pdftotext -layout) -> consensus CSV
  --compare estimates (04_quant/<pol>/nist_srm1950_estimates.tsv) + annotations -> comparison

Our estimates are summed per sum composition (one feature per unique lipid) and converted
from ng/mL to nmol/mL with the neutral mass from the feature m/z and the class's usual
positive/negative adduct ([M+NH4]+ for TG/DG/CE/MG, else [M+H]+; [M-H]- or [M+HCOO]- for
PC/SM/LPC in negative mode). A wrong adduct changes the mass by < 3 %.
"""
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import provenance  # noqa: E402

ROW = re.compile(r"^\s+(?P<name>[A-Za-z][A-Za-z0-9\-]*(?: [A-Za-z]?-?[OPd]?-?\d+:\d+[A-Za-z0-9:/_\-;()]*)?)\s{2,}"
                 r"(?P<n>\d+)\s+nmol/mL\s+(?P<dsl>[\d,.]+)\s+(?P<dsl_u>[\d,.]+)\s+(?P<vr>[\d,.]+)\s+(?P<vr_u>[\d,.]+)"
                 r"\s+(?P<medm>[\d,.]+)\s+(?P<medm_u>[\d,.]+)\s+(?P<cod>[\d,.]+)\s*$")
NAME_MAP = {"TAG": "TG", "DAG": "DG", "MAG": "MG", "FFA": "FA", "CE": "CE", "CER": "Cer"}
NH4 = {"TG", "DG", "CE", "MG"}
FORMATE_NEG = {"PC", "SM", "LPC"}
ADDUCT_MASS = {"H": 1.007276, "NH4": 18.033823, "-H": -1.007276, "HCOO": 44.997654}


def parse(text: str) -> pd.DataFrame:
    rows, table = [], None
    for line in text.splitlines():
        m = re.search(r"Table (\d+)\.", line)
        if m:
            table = int(m.group(1))
        if table is None or not 2 <= table <= 6:
            continue
        r = ROW.match(line)
        if r:
            d = r.groupdict()
            rows.append({"nist_name": d["name"].strip(), "n_labs": int(d["n"]), "table": table,
                         **{k: float(d[k].replace(",", "")) for k in ("dsl", "dsl_u", "vr", "vr_u", "medm", "medm_u", "cod")}})
    out = pd.DataFrame(rows).drop_duplicates("nist_name")
    out["sum_name"] = out["nist_name"].map(norm_nist)
    return out


def norm_nist(n: str) -> str:
    cls, _, rest = n.partition(" ")
    rest = re.sub(r"^d(?=\d)", "", rest)      # sphingolipids: 'd34:1' -> '34:1' (our sum names carry no prefix)
    return f"{NAME_MAP.get(cls, cls)} {rest}".strip()


def sum_name(lipid_name: str, lipid_class: str) -> str | None:
    """'PC 16:0_18:1' -> 'PC 34:1'; 'PC 34:1' stays; ether/odd forms -> None."""
    rest = lipid_name.split(" ", 1)[1] if " " in lipid_name else ""
    if not rest or re.search(r"[OP]-|\+|;", rest):
        return None
    chains = re.findall(r"(d|t)?(\d+):(\d+)", rest)
    if not chains:
        return None
    c, d = sum(int(x[1]) for x in chains), sum(int(x[2]) for x in chains)
    prefix = "d" if lipid_class in ("Cer", "SM") and chains[0][0] == "d" and len(chains) > 1 else ""
    if lipid_class in ("Cer", "SM"):
        return f"{lipid_class} {c}:{d}"          # NISTIR reports SM/Cer as total carbons:double bonds
    return f"{lipid_class} {prefix}{c}:{d}"


def neutral_mass(mz: float, cls: str, polarity: str) -> float:
    if polarity == "pos":
        return mz - ADDUCT_MASS["NH4" if cls in NH4 else "H"]
    return mz - ADDUCT_MASS["HCOO" if cls in FORMATE_NEG else "-H"]


def compare(est: pd.DataFrame, ann: pd.DataFrame, cons: pd.DataFrame, polarity: str) -> pd.DataFrame:
    e = est[est["n_srm"] > 0].merge(ann[["feature_id", "mz", "lipid_id"]], on="feature_id")
    e = e.sort_values("median_ng_per_ml", ascending=False).drop_duplicates("lipid_id")
    e["sum_name"] = [sum_name(n, c) for n, c in zip(e["lipid_name"], e["lipid_class"])]
    e = e.dropna(subset=["sum_name"])
    e["nmol_per_ml"] = e["median_ng_per_ml"] / [neutral_mass(m, c, polarity) for m, c in zip(e["mz"], e["lipid_class"])]
    g = e.groupby("sum_name").agg(estimate_nmol_per_ml=("nmol_per_ml", "sum"), n_features=("feature_id", "size"),
                                  routes=("quant_route", lambda s: ";".join(sorted(set(s)))),
                                  standards=("is_compound", lambda s: ";".join(sorted(set(s)))),
                                  median_rsd_pct=("rsd_pct", "median")).reset_index()
    out = g.merge(cons[["sum_name", "nist_name", "n_labs", "medm", "medm_u", "cod"]], on="sum_name")
    out["ratio"] = out["estimate_nmol_per_ml"] / out["medm"]
    out["log10_ratio"] = np.log10(out["ratio"])
    out["within_2fold"] = out["ratio"].between(0.5, 2)
    return out.sort_values("medm", ascending=False)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--parse", help="NISTIR 8185 text file")
    ap.add_argument("--consensus", required=True, help="consensus CSV (written with --parse, read otherwise)")
    ap.add_argument("--estimates")
    ap.add_argument("--annotations", help="02_annotation/<pol>/annotations.tsv (for m/z and lipid_id)")
    ap.add_argument("--polarity", choices=["pos", "neg"])
    ap.add_argument("--out")
    a = ap.parse_args(argv)
    if a.parse:
        parse(Path(a.parse).read_text()).to_csv(a.consensus, index=False)
        return
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    cons = pd.read_csv(a.consensus)
    est = pd.read_csv(a.estimates, sep="\t")
    ann = pd.read_csv(a.annotations, sep="\t")
    c = compare(est, ann, cons, a.polarity)
    c.to_csv(out / "nist_comparison.tsv", sep="\t", index=False)
    from scipy.stats import spearmanr
    rho = spearmanr(np.log10(c["estimate_nmol_per_ml"]), np.log10(c["medm"])) if len(c) > 2 else None
    summ = {"n_compared": len(c), "median_ratio": float(c["ratio"].median()) if len(c) else np.nan,
            "pct_within_2fold": float(c["within_2fold"].mean() * 100) if len(c) else np.nan,
            "spearman_rho_log": float(rho.statistic) if rho is not None else np.nan}
    for r in ("Q1", "Q3"):
        s = c[c["routes"] == r]
        summ[f"n_{r}"] = len(s)
        summ[f"median_ratio_{r}"] = float(s["ratio"].median()) if len(s) else np.nan
    pd.Series(summ).to_csv(out / "nist_summary.tsv", sep="\t", header=False)
    provenance.write(out, {"consensus": a.consensus, "estimates": a.estimates, "annotations": a.annotations},
                     {"polarity": a.polarity, "consensus_value": "MEDM location, nmol/mL, >= 5 laboratories",
                      "plasma_ml": 0.05, "aggregation": "sum per sum composition, one feature per unique lipid"})


if __name__ == "__main__":
    main()
