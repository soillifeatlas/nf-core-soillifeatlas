"""Panels: native workflow schematic with the Atlas layers, and the rendered report."""
from __future__ import annotations

import os
import shutil
import subprocess

from common import FIG_DIR, input_path, set_panel, show_image

# Text source of the native schematic below (stage names, tools, Atlas layers).
SCHEM_TEXT_SOURCE = "figures/schematics/source/Fig2_workflow_2026-11-24.png"

REPORT_HTML = "09_report/report_pos.html"
REPORT_PNG = "figures/rendered/report_pos.png"
REPORT_VIEW_CSS = (560, 300)  # narrow view: title and the summary tiles only
CHROME = "/Applications/Google Chrome.app/Contents/MacOS/Google Chrome"


def render_report() -> None:
    """Re-render the report screenshot when Chrome is available; otherwise keep the cached one."""
    html = input_path(REPORT_HTML)
    out = FIG_DIR / "rendered" / "report_pos.png"
    out.parent.mkdir(exist_ok=True)
    chrome = CHROME if os.path.exists(CHROME) else shutil.which("google-chrome")
    if chrome:
        w, h = REPORT_VIEW_CSS
        subprocess.run([chrome, "--headless=new", "--disable-gpu", "--hide-scrollbars",
                        "--force-device-scale-factor=2", f"--window-size={w},{h}",
                        f"--screenshot={out}", html.as_uri()],
                       check=True, capture_output=True)
    elif not out.exists():
        raise FileNotFoundError("no Chrome to render the report and no cached screenshot")


def report(ax):
    set_panel("report")
    render_report()
    show_image(ax, input_path(REPORT_PNG))
    ax.set_axis_on()
    ax.set_xticks([])
    ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(True)
        s.set_color("#999999")
        s.set_linewidth(0.5)


# --------------------------------------------------------------------------- native schematic
# Text of the stages, tools and Atlas layers, taken over from the 2026-11-24 draft schematic
# (Fig2_workflow_2026-11-24.png, panels A and C). Schematic labels, not data.
STAGES = [
    # (number, title, what it does, tools, output, in Nextflow?)
    (1, "Convert", "raw → indexed mzML", "ThermoRawFileParser", None, True),
    (2, "Detect", "features per polarity;\npeak-picked and\ngap-filled kept apart", "MZmine 4.9",
     None, True),
    (3, "Check quality", "pooled QCs, internal\nstandards, blanks,\ncalibration, NIST", "workflow code",
     "QC report", True),
    (4, "Annotate", "lipid names with\nconfidence levels", "LipidSearch, MS-DIAL,\nGNPS, SIRIUS, "
     "MZmine", "annotated features", True),
    (5, "Amounts", "routes Q1–Q4,\nper g dry soil", "internal standards", "amounts per g soil",
     False),
    (6, "Report", "one-page report\nand result tables", "workflow code", "report + tables",
     False),
    (7, "Interpret", "spectral match to\nreference organisms", "Atlas matching", "Atlas matches",
     False),
]
# Atlas layers under the stage they feed: (stage, text, used by default?)
ATLAS_LINKS = [(4, "identity\nreference spectra and names", False),
               (5, "response\nrelative ionisation by class", True),
               (7, "identity + occurrence\nwhere each lipid is found", True)]


def workflow_native(ax):
    """Seven stages left to right; the Atlas as a band underneath feeding stages 4, 5 and 7."""
    from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

    from common import ATLAS, SOIL
    set_panel("schematic_native")
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 40)
    ax.set_axis_off()
    x0, w, gap = 10.5, 11.3, 1.4
    top, h = 39, 22
    centres = {}
    # inputs
    ax.add_patch(FancyBboxPatch((0.3, top - h), 8.6, h, boxstyle="round,pad=0,rounding_size=1.2",
                                facecolor="#F2F2F2", edgecolor="#BBBBBB", lw=0.6))
    ax.text(4.6, top - 2.5, "Inputs", ha="center", va="top", fontsize=6.3, fontweight="bold")
    ax.text(4.6, top - 7, "raw files\n\nsample\nsheet\n\nAtlas\nrelease", ha="center", va="top",
            fontsize=5.5, linespacing=1.1)
    for i, (n, title, what, tools, out, nf) in enumerate(STAGES):
        x = x0 + i * (w + gap)
        centres[n] = x + w / 2
        ax.add_patch(FancyBboxPatch((x, top - h), w, h, boxstyle="round,pad=0,rounding_size=1.2",
                                    facecolor="white", edgecolor=SOIL, lw=0.9,
                                    linestyle="-" if nf else (0, (3, 2))))
        ax.add_patch(FancyBboxPatch((x, top - 5.2), w, 5.2,
                                    boxstyle="round,pad=0,rounding_size=1.2",
                                    facecolor=SOIL if nf else "#8A96A3", edgecolor="none"))
        ax.text(x + w / 2, top - 2.6, f"{n}  {title}", ha="center", va="center", fontsize=6,
                fontweight="bold", color="white")
        ax.text(x + w / 2, top - 7.2, what, ha="center", va="top", fontsize=5.3,
                linespacing=1.15)
        ax.text(x + w / 2, top - h + 1.2, tools, ha="center", va="bottom", fontsize=4.8,
                color="#666666", style="italic", linespacing=1.1)
        if i:
            ax.annotate("", xy=(x - 0.1, top - h / 2), xytext=(x - gap + 0.1, top - h / 2),
                        arrowprops=dict(arrowstyle="-|>", lw=0.6, color="#777777",
                                        mutation_scale=6))
        if out:
            ax.text(x + w / 2, top - h + 6.6, f"→ {out}", ha="center", va="bottom",
                    fontsize=5.2, color=SOIL, fontweight="bold")
    ax.annotate("", xy=(x0 - 0.1, top - h / 2), xytext=(9.0, top - h / 2),
                arrowprops=dict(arrowstyle="-|>", lw=0.6, color="#777777", mutation_scale=6))
    # Atlas band: each layer sits under the stage it feeds
    band_y, band_h = 0.8, 7.5
    ax.add_patch(FancyBboxPatch((x0, band_y), 7 * w + 6 * gap, band_h,
                                boxstyle="round,pad=0,rounding_size=1.2", facecolor=ATLAS,
                                alpha=0.12, edgecolor=ATLAS, lw=0.8))
    ax.text(x0 - 1.0, band_y + band_h / 2, "Soil Life\nAtlas", ha="right", va="center",
            fontsize=6.5, fontweight="bold", color=ATLAS)
    for stage, text, default in ATLAS_LINKS:
        cx = centres[stage]
        ax.text(cx, band_y + band_h / 2, text, ha="center", va="center", fontsize=5,
                color="#0B5F57", linespacing=1.15,
                bbox=dict(boxstyle="round,pad=0.3", facecolor="white", edgecolor=ATLAS, lw=0.6))
        ax.add_patch(FancyArrowPatch((cx, band_y + band_h + 0.2), (cx, top - h - 0.3),
                                     arrowstyle="-|>", mutation_scale=7, lw=0.9, color=ATLAS,
                                     linestyle="-" if default else (0, (2, 1.5))))
    ax.text(x0, band_y - 0.4, "teal arrow solid: used by default · dashed: when available    "
            "grey dashed box: run as a script, not yet a Nextflow process", ha="left", va="top",
            fontsize=5, color="#666666")
