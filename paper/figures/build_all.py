#!/usr/bin/env python
"""Rebuild every manuscript figure, STATUS.md and figure_provenance.yaml.

    python paper/figures/build_all.py

Pass-2 panels (3D, 5C, 5E-Q3, 5F, 7B-E) stay PENDING until the switches in gates.yaml are
set. Outputs go to paper/figures/output/ (or $FIG_OUT_DIR).
"""
from __future__ import annotations

import datetime as dt
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import yaml  # noqa: E402

import common  # noqa: E402
import fig2  # noqa: E402
import fig3  # noqa: E402
import fig4  # noqa: E402
import fig5  # noqa: E402
import fig6  # noqa: E402
import fig7  # noqa: E402

MODULES = [fig2, fig3, fig4, fig5, fig6, fig7]
UNCHANGED = [
    {"figure": 1, "panel": "A–D", "status": "UNCHANGED",
     "note": "not rebuilt in this pass (issue #7)", "sources": []},
    {"figure": "GA", "panel": "–", "status": "UNCHANGED",
     "note": "graphical abstract not rebuilt (issue #7)", "sources": []},
]


def _git(*args) -> str:
    return subprocess.run(["git", *args], cwd=common.REPO, capture_output=True,
                          text=True).stdout.strip()


def _export_sha(rel: str) -> str | None:
    """export_git_sha from the provenance.yaml next to a data file, if there is one."""
    folder = (common.PAPER / rel).parent
    cands = sorted(folder.glob("provenance*.yaml"))
    for c in cands:
        try:
            y = yaml.safe_load(c.read_text())
        except Exception:
            continue
        sha = (y or {}).get("export_git_sha")
        if sha:
            return str(sha)[:7]
    return None


def main() -> int:
    common.reset_reads()
    rows, outputs = [], {}
    for m in MODULES:
        paths, r = m.build()
        rows += r
        outputs[r[0]["figure"]] = [str(p.relative_to(common.PAPER)) if common.PAPER in p.parents
                                   else str(p) for p in paths]
        for p in paths:
            print("wrote", p)
    rows += UNCHANGED
    rows.sort(key=lambda r: (str(r["figure"]) if r["figure"] == "GA" else f"{r['figure']:02d}",
                             r["panel"]))

    head = _git("rev-parse", "--short", "HEAD") or "unknown"
    dirty = bool(_git("status", "--porcelain", "--", "paper/figures"))
    g = common.gates()
    src_info = {}
    for r in rows:
        for s in r["sources"]:
            if s not in src_info:
                src_info[s] = {"git_sha": common.git_sha_of(s), "export_sha": _export_sha(s),
                               "sha256": common.sha256(s)}

    # ---------------------------------------------------------------- STATUS.md
    lines = [f"# Figure status ({common.DATE})", "",
             f"Built by `python paper/figures/build_all.py` on branch "
             f"`{_git('rev-parse', '--abbrev-ref', 'HEAD')}` at `{head}`"
             f"{' (paper/figures has uncommitted changes)' if dirty else ''}.",
             "Every plotted number is read from a file under `paper/` by the figure script; "
             "nothing is typed by hand.", "",
             "Status: **REAL** = drawn from recorded data · **SCHEMATIC** = cropped from the "
             "2026-11-24 draft, no data · **PENDING** = grey box, waits for the named item · "
             "**UNCHANGED** = not rebuilt in this pass.", "",
             "Data SHA = last commit that changed the file in this repository; "
             "'export' = `export_git_sha` in the folder's provenance.yaml.", "",
             "Gates (gates.yaml): " + ", ".join(f"`{k}: {v}`" for k, v in g.items()), "",
             "| Fig | Panel | Status | Source file(s) under paper/ | Data SHA | Note |",
             "|---|---|---|---|---|---|"]
    for r in rows:
        srcs = "<br>".join(f"`{s}`" for s in r["sources"]) or "–"
        shas = "<br>".join(
            src_info[s]["git_sha"] + (f" (export {src_info[s]['export_sha']})"
                                      if src_info[s]["export_sha"] else "")
            for s in r["sources"]) or "–"
        lines.append(f"| {r['figure']} | {r['panel']} | **{r['status']}** | {srcs} | {shas} | "
                     f"{r['note']} |")
    lines += ["", "Outputs:", ""]
    for f, ps in outputs.items():
        lines.append(f"- Fig {f}: " + ", ".join(f"`{p}`" for p in ps))
    (common.FIG_DIR / "STATUS.md").write_text("\n".join(lines) + "\n")

    # ---------------------------------------------------------------- provenance
    prov = {
        "built": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),
        "command": "python paper/figures/build_all.py",
        "repo_head": head, "paper_figures_dirty": dirty, "gates": g,
        "schematic_sources": {
            "figures/schematics/source/Fig2_workflow_2026-11-24.png": {
                "from": "Google Drive folder 1bR_DKLZjM1M0crbp7-Ffwo6EUGxFjwKn, file id "
                        "1NJwtkjVm8F2GWcJSrjEDEgvCo-hli6az",
                "crops": {"2A": list(fig2.CROP_2A), "2C": list(fig2.CROP_2C)}},
            "figures/schematics/source/Fig5_quantification_2026-11-24.png": {
                "from": "Google Drive folder 1bR_DKLZjM1M0crbp7-Ffwo6EUGxFjwKn, file id "
                        "18CBhsQyLNANSEUddCtcQCL5bdV8XogRk",
                "crops": {"5D": list(fig5.CROP_5D)}},
        },
        "figures": {},
    }
    for r in rows:
        f = prov["figures"].setdefault(str(r["figure"]), {"outputs": outputs.get(r["figure"], []),
                                                          "panels": {}})
        f["panels"][r["panel"]] = {
            "status": r["status"], "note": r["note"],
            "script": f"paper/figures/fig{r['figure']}.py" if r["figure"] in outputs else None,
            "sources": {s: src_info[s] for s in r["sources"]},
        }
    (common.FIG_DIR / "figure_provenance.yaml").write_text(
        yaml.safe_dump(prov, sort_keys=False, allow_unicode=True, width=100))
    print("wrote", common.FIG_DIR / "STATUS.md")
    print("wrote", common.FIG_DIR / "figure_provenance.yaml")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
