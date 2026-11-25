"""Panels: workflow schematic (crop), tools/Atlas-layer matrix (crop), rendered report."""
from __future__ import annotations

import os
import shutil
import subprocess

from matplotlib.patches import FancyBboxPatch, Rectangle

from common import FIG_DIR, input_path, set_panel, show_image

SCHEM = "figures/schematics/source/Fig2_workflow_2026-11-24.png"
CROP_2A = (102, 240, 3866, 944)    # px in the 4000x3160 draft
CROP_2C = (2183, 1208, 3885, 2100)
STAGE_BOXES = [(2358, 2792), (2834, 3268), (3310, 3745)]  # stages 5, 6, 7 inside CROP_2A
STAGE_Y = (18, 448)

REPORT_HTML = "09_report/report_pos.html"
REPORT_PNG = "figures/rendered/report_pos.png"
REPORT_VIEW_CSS = (1100, 470)  # top of the page: title, summary tiles, first quality checks
CHROME = "/Applications/Google Chrome.app/Contents/MacOS/Google Chrome"


def stages(ax):
    """Seven-stage schematic; stages 5-7 dashed (scripted, not yet Nextflow processes)."""
    set_panel("schematic_stages")
    show_image(ax, input_path(SCHEM), CROP_2A)
    y0, y1 = STAGE_Y
    for x0, x1 in STAGE_BOXES:
        ax.add_patch(Rectangle((x0, y0), x1 - x0, y1 - y0, facecolor="white", alpha=0.55,
                               edgecolor="none", zorder=3))
        ax.add_patch(FancyBboxPatch((x0 + 6, y0 + 6), x1 - x0 - 12, y1 - y0 - 12,
                                    boxstyle="round,pad=0,rounding_size=30", fill=False,
                                    edgecolor="#333333", linewidth=0.9,
                                    linestyle=(0, (4, 2.5)), zorder=4))


def tools_matrix(ax):
    set_panel("schematic_tools")
    show_image(ax, input_path(SCHEM), CROP_2C)


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
