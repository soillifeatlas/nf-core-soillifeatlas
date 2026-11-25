"""Figure 2 - workflow. A, C schematic crops; B real sample sheet; D rendered stage-6 report."""
from __future__ import annotations

import shutil
import subprocess

import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Rectangle

from common import (DOUBLE, FIG_DIR, MM, ROLE, ROLE_LABEL, Panels, fig_letter, input_path,
                    read_table, save, set_panel, show_image)

SCHEM = "figures/schematics/source/Fig2_workflow_2026-11-24.png"
CROP_2A = (102, 240, 3866, 944)    # px in the 4000x3160 draft
CROP_2C = (2183, 1208, 3885, 2100)
# stage boxes 5, 6, 7 inside the 2A crop (px): x0, x1; shared y0, y1
STAGE_BOXES = [(2358, 2792), (2834, 3268), (3310, 3745)]
STAGE_Y = (18, 448)

REPORT_HTML = "09_report/report_pos.html"
REPORT_PNG = "figures/rendered/report_pos.png"
CHROME = "/Applications/Google Chrome.app/Contents/MacOS/Google Chrome"


def render_report() -> None:
    """Re-render the report screenshot when Chrome is available; otherwise keep the cached one."""
    html = input_path(REPORT_HTML)
    out = FIG_DIR / "rendered" / "report_pos.png"
    out.parent.mkdir(exist_ok=True)
    chrome = CHROME if shutil.os.path.exists(CHROME) else shutil.which("google-chrome")
    if chrome:
        subprocess.run([chrome, "--headless=new", "--disable-gpu", "--hide-scrollbars",
                        "--force-device-scale-factor=2", "--window-size=1100,745",
                        f"--screenshot={out}", html.as_uri()],
                       check=True, capture_output=True)
    elif not out.exists():
        raise FileNotFoundError("no Chrome to render the report and no cached screenshot")


def build():
    P = Panels(2)
    fig = plt.figure(figsize=(DOUBLE, 225 * MM))

    # ---- A: schematic crop, stages 5-7 marked as scripts
    set_panel("2A")
    axA = fig.add_axes([0.0, 0.815, 1.0, 0.155])
    show_image(axA, input_path(SCHEM), CROP_2A)
    y0, y1 = STAGE_Y
    for x0, x1 in STAGE_BOXES:
        axA.add_patch(Rectangle((x0, y0), x1 - x0, y1 - y0, facecolor="white", alpha=0.55,
                                edgecolor="none", zorder=3))
        axA.add_patch(FancyBboxPatch((x0 + 6, y0 + 6), x1 - x0 - 12, y1 - y0 - 12,
                                     boxstyle="round,pad=0,rounding_size=30", fill=False,
                                     edgecolor="#333333", linewidth=0.9,
                                     linestyle=(0, (4, 2.5)), zorder=4))
        axA.text((x0 + x1) / 2, y1 - 60, "script, not yet\nin Nextflow", ha="center",
                 va="center", fontsize=5.5, fontstyle="italic", color="#222222", zorder=5)
    fig_letter(fig, "A", 0.005, 0.995,
               "Workflow stages (1–4 run in Nextflow; 5–7 run as library code and scripts)")
    P.add("A", "SCHEMATIC", "crop of Fig2_workflow_2026-11-24.png; stages 5–7 overlaid dashed "
          "('script, not yet in Nextflow')")

    # ---- B: real sample sheet rows + injection strip
    set_panel("2B")
    ss = read_table("09_report/samplesheet_climgrass_pos.csv")
    inj = read_table("01_qc/pos/injections.tsv")
    ss = ss.sort_values("injection_order").reset_index(drop=True)
    cal = ss[ss.role == "standard_cal"].head(2)
    rest = ss[ss.role != "standard_cal"].head(9)
    shown = [cal, None, rest]
    axB = fig.add_axes([0.03, 0.50, 0.52, 0.27])
    axB.set_axis_off()
    fig_letter(fig, "B", 0.005, 0.79, "Sample sheet and injection sequence (ClimGrass, POS)")
    cols = [("file", "raw file", 0.0, "left"), ("role", "role", 0.37, "left"),
            ("injection_order", "order", 0.62, "right"), ("group", "treatment", 0.65, "left"),
            ("dry_mass_g", "dry mass (g)", 0.92, "right"), ("spike_ng", "IS (ng)", 1.0, "right")]
    y = 1.0
    dy = 0.062
    for key, head, x, ha in cols:
        axB.text(x, y, head, fontweight="bold", ha=ha, va="top", transform=axB.transAxes,
                 fontsize=6)
    axB.plot([0, 1], [y - 0.045, y - 0.045], color="#888888", lw=0.5, transform=axB.transAxes)
    y -= dy + 0.01
    for block in shown:
        if block is None:
            axB.text(0.0, y, "…", va="top", transform=axB.transAxes, color="#666666")
            y -= dy
            continue
        for _, r in block.iterrows():
            for key, _, x, ha in cols:
                v = r[key]
                s = "–" if v != v or v == "" else (ROLE_LABEL.get(v, v) if key == "role" else
                                                   f"{v:.3f}" if key == "dry_mass_g" else
                                                   f"{v:g}" if key == "spike_ng" else
                                                   str(v).replace("_", " + ") if key == "group"
                                                   else str(v))
                axB.text(x, y, s, ha=ha, va="top", transform=axB.transAxes, fontsize=6)
            y -= dy
    axB.text(0.0, y, f"…  {len(ss)} injections in total, in run order",
             va="top", transform=axB.transAxes, color="#666666", fontsize=6)

    axS = fig.add_axes([0.03, 0.455, 0.52, 0.025])
    blocks = dict(zip(inj.injection_order, inj.block))
    for _, r in ss.iterrows():
        o = r.injection_order
        axS.add_patch(Rectangle((o - 0.45, 0), 0.9, 1, facecolor=ROLE[r.role],
                                edgecolor="#333333", linewidth=0.3))
    b2 = min(o for o, b in blocks.items() if b == 2)
    axS.axvline(b2 - 0.5, color="black", lw=0.8)
    axS.set_xlim(0.4, len(ss) + 0.6)
    axS.set_ylim(0, 1)
    axS.set_yticks([])
    axS.set_xticks([1, b2, len(ss)])
    axS.spines[["left"]].set_visible(False)
    axS.set_xlabel("injection order", labelpad=1)
    axS.text(b2 / 2, 1.25, "block 1: calibration", ha="center", fontsize=6)
    axS.text((b2 + len(ss)) / 2, 1.25, "block 2: soils, QCs, SRM", ha="center", fontsize=6)
    handles = [Rectangle((0, 0), 1, 1, facecolor=ROLE[k], edgecolor="#333333", linewidth=0.3)
               for k in ROLE]
    fig.legend(handles, [ROLE_LABEL[k] for k in ROLE], loc="upper left",
               bbox_to_anchor=(0.03, 0.425), ncol=3, handlelength=1.0, columnspacing=1.0)
    P.add("B", "REAL", "first two calibration rows and the first nine soil-block rows shown")

    # ---- C: schematic crop
    set_panel("2C")
    axC = fig.add_axes([0.58, 0.47, 0.41, 0.30])
    show_image(axC, input_path(SCHEM), CROP_2C)
    fig_letter(fig, "C", 0.565, 0.79, "Where tools and Atlas layers enter")
    P.add("C", "SCHEMATIC", "crop of Fig2_workflow_2026-11-24.png, unchanged (#10)")

    # ---- D: rendered report
    set_panel("2D")
    render_report()
    axD = fig.add_axes([0.03, 0.0, 0.94, 0.345])
    show_image(axD, input_path(REPORT_PNG))
    axD.set_axis_on()
    axD.set_xticks([])
    axD.set_yticks([])
    for s in axD.spines.values():
        s.set_visible(True)
        s.set_color("#999999")
        s.set_linewidth(0.5)
    fig_letter(fig, "D", 0.005, 0.365, "Stage-6 report (ClimGrass, POS), top of page")
    P.add("D", "REAL", "headless-Chrome render of report_pos.html (top section); the "
          "'Lipids with an amount' tile is stage-5 output")
    return save(fig, 2, "workflow"), P.rows
