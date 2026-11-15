#!/usr/bin/env python3
"""Stage 3: quality-check report for one polarity of a processed run.

Inputs: the MZmine full feature table, the sample sheet, the internal-standard
target table and (optionally) the mzML files, whose headers give the injection
times. Outputs, all prefixed with --prefix:

  <prefix>.qc_checks.tsv               pass / warn / fail per check
  <prefix>.qc_internal_standards.tsv   one row per EquiSPLASH component
  <prefix>.qc_is_per_injection.tsv     standard area, RT and ppm per injection
  <prefix>.qc_injections.tsv           features and total area per injection
  <prefix>.qc_features.tsv             pooled-QC RSD and blank ratio per feature
  <prefix>.qc_summary.json             thresholds and overall status
  <prefix>.qc_report.html              the same, readable, with plots
"""
from __future__ import annotations

import argparse
import base64
import html
import io
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from soillifeatlas.process_qc import (
    Thresholds,
    file_id,
    read_mzml_start_time,
    read_samplesheet,
    read_targets,
    run_qc,
)

ROLE_COLOURS = {
    "sample": "#2f6b4f",
    "pooled_qc": "#e0a85e",
    "blank": "#9aa89f",
    "extraction_blank": "#6a786f",
    "nist": "#7b4fa0",
    "standard_cal": "#245276",
    "not_in_sheet": "#c0392b",
}


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--features", required=True, help="MZmine full feature table (CSV)")
    p.add_argument("--samplesheet", required=True)
    p.add_argument("--targets", required=True, help="Internal-standard targets (CSV)")
    p.add_argument("--polarity", required=True, choices=["pos", "neg"])
    p.add_argument("--mzml", nargs="*", default=[], help="mzML files, for injection times")
    p.add_argument("--prefix", required=True)
    p.add_argument("--rt-tolerance", type=float, default=0.1)
    p.add_argument("--is-rsd", type=float, default=20.0)
    p.add_argument("--mass-accuracy-ppm", type=float, default=5.0)
    p.add_argument("--qc-rsd", type=float, default=30.0)
    p.add_argument("--qc-pass-fraction", type=float, default=0.7)
    p.add_argument("--blank-ratio", type=float, default=0.2)
    p.add_argument("--match-ppm", type=float, default=10.0)
    p.add_argument("--rt-window", type=float, default=0.3,
                   help="Window (min) around an expected RT given in the standards table")
    p.add_argument("--is-compounds", default="",
                   help="Comma-separated subset of standards to use (default: all in the table)")
    return p.parse_args(argv)


def overall_status(checks: pd.DataFrame) -> str:
    if (checks["status"] == "fail").any():
        return "fail"
    if (checks["status"] == "warn").any():
        return "warn"
    return "pass"


# ---------------------------------------------------------------------------
# Plots (optional: the report is still written without matplotlib)
# ---------------------------------------------------------------------------
def _fig_to_img(fig) -> str:
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=110, bbox_inches="tight")
    import matplotlib.pyplot as plt

    plt.close(fig)
    return '<img alt="" src="data:image/png;base64,' + base64.b64encode(buf.getvalue()).decode() + '">'


def make_plots(res) -> dict[str, str]:
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        return {}
    plots: dict[str, str] = {}
    inj = res.injections
    order = inj.set_index("file_id")["injection_order"]

    fig, ax = plt.subplots(figsize=(9, 3.2))
    for role, grp in inj.groupby("role"):
        ax.scatter(grp["injection_order"], grp["features_detected"], s=28,
                   color=ROLE_COLOURS.get(role, "#333"), label=role)
    ax.set_xlabel("Injection order")
    ax.set_ylabel("Features detected")
    ax.legend(fontsize=8, ncol=3, frameon=False)
    ax.grid(alpha=0.25)
    plots["injections"] = _fig_to_img(fig)

    det = res.is_per_injection
    if det is not None and not det.empty:
        det = det.copy()
        det["injection_order"] = det["file_id"].map(order)
        det = det.sort_values("injection_order")
        fig, axes = plt.subplots(1, 3, figsize=(13, 3.4))
        for compound, grp in det.groupby("compound"):
            g = grp[grp["role"].isin(["sample", "pooled_qc"])]
            if g["area"].notna().sum() < 2:
                continue
            axes[0].plot(g["injection_order"], g["rt"] - g["rt"].median(), marker="o", ms=3, lw=1, label=compound)
            axes[1].plot(g["injection_order"], g["ppm_error"], marker="o", ms=3, lw=1)
            axes[2].plot(g["injection_order"], g["area"] / g["area"].median(), marker="o", ms=3, lw=1)
        for ax, title in zip(axes, ["RT minus median (min)", "Mass error (ppm)", "Area / median area"]):
            ax.set_title(title, fontsize=10)
            ax.set_xlabel("Injection order")
            ax.grid(alpha=0.25)
        axes[2].set_yscale("log")
        fig.legend(*axes[0].get_legend_handles_labels(), loc="lower center", ncol=5, fontsize=7, frameon=False,
                   bbox_to_anchor=(0.5, -0.12))
        plots["internal_standards"] = _fig_to_img(fig)

    feats = res.features
    if feats is not None and "qc_rsd_pct" in feats:
        vals = feats.loc[feats["qc_detected"] >= 2, "qc_rsd_pct"].dropna().clip(upper=200)
        if len(vals):
            fig, ax = plt.subplots(figsize=(6, 3))
            ax.hist(vals, bins=40, color="#e0a85e")
            ax.axvline(30, color="#1d3a2c", ls="--", lw=1)
            ax.set_xlabel("Pooled-QC RSD (%) per feature (capped at 200)")
            ax.set_ylabel("Features")
            plots["qc_rsd"] = _fig_to_img(fig)
    return plots


# ---------------------------------------------------------------------------
# HTML
# ---------------------------------------------------------------------------
STATUS_STYLE = {
    "pass": "background:#e2efe7;color:#1d3a2c",
    "warn": "background:#f6ead8;color:#7a4d14",
    "fail": "background:#f6dcd8;color:#8a1c1c",
    "not_run": "background:#eceeed;color:#4f5d54",
}


def _table(df: pd.DataFrame, float_fmt="{:.3g}") -> str:
    def fmt(v):
        if isinstance(v, float):
            return "" if np.isnan(v) else float_fmt.format(v)
        return html.escape(str(v))

    head = "".join(f"<th>{html.escape(str(c))}</th>" for c in df.columns)
    body = []
    for _, r in df.iterrows():
        cells = []
        for c in df.columns:
            style = f' style="{STATUS_STYLE[r[c]]}"' if c == "status" and r[c] in STATUS_STYLE else ""
            cells.append(f"<td{style}>{fmt(r[c])}</td>")
        body.append("<tr>" + "".join(cells) + "</tr>")
    return f'<div class="tw"><table><thead><tr>{head}</tr></thead><tbody>{"".join(body)}</tbody></table></div>'


def write_html(path: Path, res, checks: pd.DataFrame, summary: dict, plots: dict[str, str]) -> None:
    status = summary["overall_status"]
    is_cols = [c for c in ["compound", "lipid_class", "matched", "adduct", "priority", "selected_by", "median_rt",
                           "max_rt_deviation_min", "median_ppm_error", "max_abs_ppm_error", "area_rsd_pct",
                           "detected_in", "comparable_injections"] if c in res.internal_standards.columns]
    doc = f"""<!doctype html><html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>QC report · {html.escape(summary['prefix'])}</title>
<style>
body{{font-family:Inter,system-ui,sans-serif;margin:0;background:#f4f6f2;color:#18231d;line-height:1.5}}
main{{max-width:1100px;margin:0 auto;padding:24px 16px 48px}}
h1{{font-size:26px;margin:0 0 4px}} h2{{font-size:19px;margin:32px 0 8px}}
.badge{{display:inline-block;padding:4px 12px;border-radius:999px;font-weight:600;{STATUS_STYLE.get(status, '')}}}
.tw{{overflow-x:auto;background:#fff;border:1px solid #d5ded8;border-radius:8px}}
table{{border-collapse:collapse;width:100%;font-size:13px}}
th,td{{text-align:left;padding:6px 10px;border-bottom:1px solid #e1e7e3;vertical-align:top}}
th{{font-size:11px;text-transform:uppercase;letter-spacing:.05em;color:#4f5d54}}
img{{max-width:100%;background:#fff;border:1px solid #d5ded8;border-radius:8px}}
.meta{{color:#4f5d54;font-size:13px}}
</style></head><body><main>
<h1>Quality check · {html.escape(summary['prefix'])}</h1>
<p class="meta">Polarity {summary['polarity']} · {summary['n_injections']} injections · {summary['n_features']} features ·
injection order from {html.escape(summary['order_source'])}</p>
<p>Overall: <span class="badge">{status}</span></p>
<h2>Checks</h2>
{_table(checks)}
<h2>Internal standards</h2>
<p class="meta">Retention time, mass error and area RSD are computed over samples and pooled QCs only.</p>
{_table(res.internal_standards[is_cols])}
{plots.get('internal_standards', '')}
<h2>Injections</h2>
{plots.get('injections', '')}
{_table(res.injections)}
<h2>Pooled-QC reproducibility</h2>
{plots.get('qc_rsd', '<p class="meta">Needs at least two pooled QCs.</p>')}
<p class="meta">Thresholds: {html.escape(json.dumps(summary['thresholds']))}</p>
</main></body></html>"""
    path.write_text(doc)


def main(argv=None) -> int:
    a = parse_args(argv)
    full = pd.read_csv(a.features, low_memory=False)
    sheet = read_samplesheet(a.samplesheet)
    targets = read_targets(a.targets, [c for c in a.is_compounds.split(",") if c.strip()])
    thresholds = Thresholds(
        rt_tolerance_min=a.rt_tolerance,
        is_area_rsd_pct=a.is_rsd,
        mass_accuracy_ppm=a.mass_accuracy_ppm,
        qc_feature_rsd_pct=a.qc_rsd,
        qc_feature_pass_fraction=a.qc_pass_fraction,
        blank_ratio=a.blank_ratio,
        match_ppm=a.match_ppm,
        match_rt_window_min=a.rt_window,
    )
    start_times = {file_id(p): read_mzml_start_time(p) for p in a.mzml}
    res = run_qc(full, sheet, targets, a.polarity, start_times, thresholds)
    checks = res.checks_frame()

    out = Path(".")
    checks.to_csv(out / f"{a.prefix}.qc_checks.tsv", sep="\t", index=False)
    res.internal_standards.to_csv(out / f"{a.prefix}.qc_internal_standards.tsv", sep="\t", index=False)
    res.is_per_injection.to_csv(out / f"{a.prefix}.qc_is_per_injection.tsv", sep="\t", index=False)
    res.injections.to_csv(out / f"{a.prefix}.qc_injections.tsv", sep="\t", index=False)
    res.features.to_csv(out / f"{a.prefix}.qc_features.tsv", sep="\t", index_label="feature_id")

    summary = {
        "prefix": a.prefix,
        "polarity": a.polarity,
        "overall_status": overall_status(checks),
        "n_injections": int(len(res.injections)),
        "n_features": int(len(res.features)),
        "order_source": str(res.injections["order_source"].iloc[0]) if len(res.injections) else "-",
        "thresholds": thresholds.__dict__,
        "checks": checks.to_dict(orient="records"),
    }
    (out / f"{a.prefix}.qc_summary.json").write_text(json.dumps(summary, indent=2, default=str))
    write_html(out / f"{a.prefix}.qc_report.html", res, checks, summary, make_plots(res))
    print(f"QC {a.prefix}: {summary['overall_status']} "
          f"({(checks['status'] == 'fail').sum()} fail, {(checks['status'] == 'warn').sum()} warn)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
