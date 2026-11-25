"""Shared helpers for the manuscript figures.

Rules enforced here:
  * every data file is read through `read_table()` / `read_yaml()`, which refuse paths
    outside paper/ and record what each table or panel read (for STATUS.md and provenance);
  * numbers are made only in analysis.py; plots.py and captions.py read its outputs;
  * style: Elsevier widths (90 / 190 mm), Arial 7 pt, bold 9 pt panel letters,
    Okabe-Ito colour-blind-safe palette, POS/NEG colours fixed for all figures.
"""
from __future__ import annotations

import hashlib
import os
import subprocess
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import pandas as pd  # noqa: E402
import yaml  # noqa: E402

FIG_DIR = Path(__file__).resolve().parent
PAPER = FIG_DIR.parent
REPO = PAPER.parent
OUT_DIR = Path(os.environ.get("FIG_OUT_DIR", FIG_DIR / "output"))
DATE = "2026-11-25"  # file-name stamp only

MM = 1 / 25.4
SINGLE = 90 * MM
DOUBLE = 190 * MM

# Okabe-Ito
OI = {
    "black": "#000000", "orange": "#E69F00", "sky": "#56B4E9", "green": "#009E73",
    "yellow": "#F0E442", "blue": "#0072B2", "vermillion": "#D55E00", "purple": "#CC79A7",
    "grey": "#999999",
}
POL = {"pos": OI["blue"], "neg": OI["vermillion"]}
POL_LABEL = {"pos": "POS", "neg": "NEG"}
ROLE = {
    "standard_cal": OI["orange"], "pooled_qc": OI["green"], "sample": OI["sky"],
    "nist": OI["purple"], "blank": "#FFFFFF", "extraction_blank": OI["yellow"],
}
ROLE_LABEL = {
    "standard_cal": "calibration standard", "pooled_qc": "pooled QC", "sample": "soil",
    "nist": "NIST SRM 1950", "blank": "solvent blank", "extraction_blank": "extraction blank",
}
CONF = {"Gold": OI["orange"], "Silver": OI["grey"], "Bronze": "#8C510A"}

plt.rcParams.update({
    "font.family": "Arial", "font.size": 7, "axes.titlesize": 7, "axes.labelsize": 7,
    "xtick.labelsize": 6.5, "ytick.labelsize": 6.5, "legend.fontsize": 6.5,
    "axes.linewidth": 0.6, "xtick.major.width": 0.6, "ytick.major.width": 0.6,
    "xtick.major.size": 2.5, "ytick.major.size": 2.5, "lines.linewidth": 0.9,
    "axes.spines.top": False, "axes.spines.right": False, "legend.frameon": False,
    "pdf.fonttype": 42, "mathtext.fontset": "custom", "mathtext.rm": "Arial",
    "mathtext.it": "Arial:italic", "mathtext.bf": "Arial:bold", "mathtext.default": "regular", "ps.fonttype": 42, "savefig.dpi": 600, "figure.dpi": 150,
})

# --------------------------------------------------------------------------- reading
_READS: dict[str, set[str]] = {}
_current_panel: list[str] = [""]


def set_panel(key: str) -> None:
    """Attribute subsequent reads to a panel key such as '3B'."""
    _current_panel[0] = key
    _READS.setdefault(key, set())


def reads() -> dict[str, list[str]]:
    return {k: sorted(v) for k, v in _READS.items()}


def reset_reads() -> None:
    _READS.clear()


def _resolve(rel: str) -> Path:
    p = (PAPER / rel).resolve()
    if PAPER.resolve() not in p.parents:
        raise PermissionError(f"figure input outside paper/: {rel}")
    if not p.exists():
        raise FileNotFoundError(f"expected input missing: paper/{rel}")
    if _current_panel[0]:
        _READS[_current_panel[0]].add(str(p.relative_to(PAPER)))
    return p


def read_table(rel: str, **kw) -> pd.DataFrame:
    p = _resolve(rel)
    sep = "," if p.suffix == ".csv" else "\t"
    return pd.read_csv(p, sep=sep, **kw)


def read_yaml(rel: str) -> dict:
    return yaml.safe_load(_resolve(rel).read_text())


def input_path(rel: str) -> Path:
    """Path of a non-table input (image, HTML) inside paper/, recorded like a read."""
    return _resolve(rel)


def gates() -> dict:
    """Release switches for Pass-2 panels (paper/figures/gates.yaml)."""
    path = Path(os.environ.get("FIG_GATES_FILE", FIG_DIR / "gates.yaml"))
    return yaml.safe_load(path.read_text())


# --------------------------------------------------------------------------- provenance
def git_sha_of(rel: str) -> str:
    """Short SHA of the last commit that changed paper/<rel> (the data version)."""
    try:
        out = subprocess.run(["git", "log", "-1", "--format=%h", "--", f"paper/{rel}"],
                             cwd=REPO, capture_output=True, text=True, check=True).stdout.strip()
        dirty = subprocess.run(["git", "status", "--porcelain", "--", f"paper/{rel}"],
                               cwd=REPO, capture_output=True, text=True).stdout.strip()
        return (out or "uncommitted") + ("+dirty" if dirty else "")
    except Exception:
        return "unknown"


def sha256(rel: str) -> str:
    return hashlib.sha256((PAPER / rel).read_bytes()).hexdigest()


# --------------------------------------------------------------------------- drawing
def letter(ax, s: str, x: float = -0.02, y: float = 1.0, title: str | None = None) -> None:
    """Bold 9 pt panel letter at the top-left of an axes (axes-fraction coordinates)."""
    ax.text(x, y, s, transform=ax.transAxes, fontsize=9, fontweight="bold",
            va="bottom", ha="right")
    if title:
        ax.text(x + 0.015, y, title, transform=ax.transAxes, fontsize=7,
                va="bottom", ha="left")


def fig_letter(fig, s: str, x: float, y: float, title: str | None = None) -> None:
    """Panel letter in figure-fraction coordinates (for image/schematic panels)."""
    fig.text(x, y, s, fontsize=9, fontweight="bold", va="top", ha="left")
    if title:
        fig.text(x + 0.018, y - 0.001, title, fontsize=7, va="top", ha="left")


def show_image(ax, path: Path, crop: tuple[int, int, int, int] | None = None):
    from PIL import Image
    im = Image.open(path).convert("RGB")
    if crop:
        im = im.crop(crop)
    ax.imshow(im, interpolation="lanczos")
    ax.set_axis_off()
    return im


def save(fig, number: int, name: str) -> list[Path]:
    OUT_DIR.mkdir(exist_ok=True)
    stem = OUT_DIR / f"Fig{number}_{name}_{DATE}"
    paths = []
    for ext in ("pdf", "png"):
        p = stem.with_suffix(f".{ext}")
        fig.savefig(p, dpi=600, metadata={"Creator": "paper/figures", "CreationDate": None}
                    if ext == "pdf" else None)
        paths.append(p)
    plt.close(fig)
    return paths


class Panels:
    """Collects one STATUS row per panel. `keys` name the set_panel() buckets it read from."""

    def __init__(self, figure: str, name: str, role: str = "main"):
        self.figure, self.name, self.role = str(figure), name, role
        self.rows: list[dict] = []

    def add(self, panel: str, status: str, note: str = "", keys: list[str] | None = None):
        r = reads()
        src = sorted({f for k in (keys or []) for f in r.get(k, [])})
        self.rows.append({"figure": self.figure, "name": self.name, "role": self.role,
                          "panel": panel, "status": status, "sources": src, "note": note})


# --------------------------------------------------------------------------- shared vocabulary
# One meaning per colour in the main figures: soil = slate, Atlas = teal, named/annotated =
# mid grey, no match = light grey; organism groups use Paul Tol's colour-blind-safe set.
SOIL = "#3B4A5A"
ATLAS = "#009988"
NAMED = "#8C8C8C"
NOMATCH = "#DADADA"
GROUP_ORDER = ["Bacteria", "Archaea", "Fungi", "Protists", "Plants", "Soil animals"]
GROUP_COLOR = {"Bacteria": "#4477AA", "Archaea": "#AA3377", "Fungi": "#CCBB44",
               "Protists": "#66CCEE", "Plants": "#117733", "Soil animals": "#EE6677"}


def atlas_groups() -> pd.DataFrame:
    """Phylum -> organism group, exported from the release by export_groups.py."""
    return read_table("figures/inputs/atlas_groups.csv")


def finding(ax, tag: str, text: str, x: float = 0.0, y: float = 1.04) -> None:
    """Panel letter followed by a one-line finding (the take-home message)."""
    ax.text(x, y, tag, transform=ax.transAxes, fontsize=9, fontweight="bold", va="bottom",
            ha="right")
    ax.text(x + 0.012, y, text, transform=ax.transAxes, fontsize=7, va="bottom", ha="left")
