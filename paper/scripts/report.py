#!/usr/bin/env python3
"""Issue #9: minimal stage-6 report (one HTML page) and the real ClimGrass sample sheet.

Reads the paper exports (01_qc, 02_annotation, 04_quant, 08_application) of one
polarity and writes report_<pol>.html plus samplesheet_climgrass_<pol>.csv with
the pipeline's sample-sheet columns (file, polarity, role, group, dry_mass_g, spike_ng).
"""
from __future__ import annotations

import argparse
import html
import re
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import provenance  # noqa: E402

CSS = """body{font:14px/1.45 system-ui,sans-serif;margin:24px auto;max-width:1000px;color:#1f2a24;background:#fff}
h1{font-size:22px;margin:0 0 4px}h2{font-size:16px;margin:22px 0 6px;border-bottom:1px solid #d6ddd8}
table{border-collapse:collapse;font-size:12.5px;margin:4px 0}td,th{padding:3px 8px;border-bottom:1px solid #eef1ef;text-align:right}
th{background:#f3f6f4}td:first-child,th:first-child{text-align:left}.pass{color:#1d7a46}.warn{color:#a86b00}.fail{color:#b3261e}
.cards{display:flex;gap:10px;flex-wrap:wrap}.card{border:1px solid #d6ddd8;border-radius:6px;padding:8px 12px;min-width:150px}
.card b{display:block;font-size:20px}.muted{color:#66736b;font-size:12px}"""


def table(df: pd.DataFrame, n: int | None = None) -> str:
    d = df.head(n) if n else df
    head = "".join(f"<th>{html.escape(str(c))}</th>" for c in d.columns)
    body = ""
    for _, r in d.iterrows():
        cells = []
        for c in d.columns:
            v = r[c]
            s = f"{v:,.3g}" if isinstance(v, float) else html.escape(str(v))
            cls = f' class="{v}"' if c == "status" and v in ("pass", "warn", "fail") else ""
            cells.append(f"<td{cls}>{s}</td>")
        body += "<tr>" + "".join(cells) + "</tr>"
    return f"<table><tr>{head}</tr>{body}</table>"


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--paper", required=True, help="paper/ directory")
    ap.add_argument("--polarity", choices=["pos", "neg"], required=True)
    ap.add_argument("--samples", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args(argv)
    P, pol, out = Path(a.paper), a.polarity, Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    f = {"checks": P / "01_qc" / pol / "checks.tsv", "rsd": P / "01_qc" / pol / "rsd_summary.tsv",
         "injections": P / "01_qc" / pol / "injections.tsv", "confidence": P / "02_annotation" / pol / "confidence.tsv",
         "routes": P / "04_quant" / pol / "route_counts.tsv", "profiles": P / "08_application" / pol / "class_profiles.tsv",
         "responsive": P / "08_application" / pol / "responsive.tsv", "per_lipid": P / "04_quant" / pol / "per_lipid.tsv",
         "samples": Path(a.samples)}
    checks = pd.read_csv(f["checks"], sep="\t")
    rsd = pd.read_csv(f["rsd"], sep="\t")
    conf = pd.read_csv(f["confidence"], sep="\t")
    routes = pd.read_csv(f["routes"], sep="\t")
    prof = pd.read_csv(f["profiles"], sep="\t")
    wide = prof.pivot(index="lipid_class", columns="treatment", values="mean")
    wide = wide.loc[wide.mean(axis=1).sort_values(ascending=False).index].reset_index().head(12)
    resp = pd.read_csv(f["responsive"], sep="\t")
    per = pd.read_csv(f["per_lipid"], sep="\t")
    top = per[per["quant_route"].isin(["Q1", "Q2", "Q3"])].sort_values("median_ng_per_g", ascending=False)
    top = top[["feature_id", "lipid_name", "confidence", "quant_route", "is_compound", "f", "f_capped",
               "median_ng_per_g"]].head(12)
    top["median_ug_per_g"] = top.pop("median_ng_per_g") / 1000

    n_pass = int((checks["status"] == "pass").sum())
    pp = rsd.set_index("mode").loc["peak_picked_only"]
    cards = [("QC checks passed", f"{n_pass} / {len(checks)}"),
             ("Pooled-QC median RSD", f"{float(pp['median_rsd_pct']):.1f} %"),
             ("Annotated features", f"{int(conf.loc[conf.confidence != 'unidentified', 'n_features'].sum()):,}"),
             ("Lipids with an amount (Q1-Q3)", f"{int(routes.loc[routes.quant_route.isin(['Q1','Q2','Q3']), 'n_unique_lipids'].sum()):,}")]
    page = f"""<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>ClimGrass run report</title><style>{CSS}</style></head><body>
<h1>ClimGrass · {'positive' if pol == 'pos' else 'negative'} mode</h1>
<div class="muted">Soil Life Atlas workflow, stages 1–5 · 12 soils, 2 × 2 design (drought × +3 °C/+300 ppm CO₂), 3 plots each</div>
<h2>At a glance</h2><div class="cards">{''.join(f'<div class="card"><span class="muted">{k}</span><b>{v}</b></div>' for k, v in cards)}</div>
<h2>Quality checks</h2>{table(checks[['check', 'threshold', 'value', 'status']])}
<h2>Annotation confidence</h2>{table(conf)}
<h2>Quantification routes</h2>{table(routes)}
<h2>Class amounts, µg g⁻¹ dry soil (mean of 3 plots)</h2>{table(wide)}
<h2>Largest lipid amounts</h2>{table(top)}
<h2>Drought-responsive lipids (two-way ANOVA, BH q &lt; 0.05)</h2>{table(resp[['lipid_name', 'lipid_class', 'confidence', 'quant_route', 'log2fc_drought_vs_ambient', 'q_drought']], 15)}
<p class="muted">Q3 amounts use class-level responses (Samrat et al. 2025, Table S10) with the response factor limited to 0.1–10; capped rows are marked. Amounts are semi-quantitative.</p>
</body></html>"""
    (out / f"report_{pol}.html").write_text(page)

    # real sample sheet rows (Fig. 2B)
    inj = pd.read_csv(f["injections"], sep="\t").sort_values("injection_order")
    meta = pd.read_csv(f["samples"])
    plot = inj["file_id"].map(lambda s: int(m.group(1)) if (m := re.match(r"OE12-3-(\d+)(POS|NEG)$", s)) else None)
    mm = meta.set_index("plot")
    sheet = pd.DataFrame({"file": inj["file_id"] + ".raw", "polarity": pol, "role": inj["role"],
                          "injection_order": inj["injection_order"],
                          "group": plot.map(mm["treatment"]), "dry_mass_g": plot.map(mm["dry_mass_g"]),
                          "spike_ng": plot.map(mm["spike_ng_per_is"])})
    sheet.to_csv(out / f"samplesheet_climgrass_{pol}.csv", index=False)
    provenance.write(out, f, {"polarity": pol})


if __name__ == "__main__":
    main()
