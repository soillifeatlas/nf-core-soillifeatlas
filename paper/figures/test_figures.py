"""Checks for the figure code.  Run:  python paper/figures/test_figures.py  (or pytest).

1. analysis.py reproduces expected_numbers.csv (every number behind a figure or caption).
2. Every figure builds with the gates as committed; gated figures are reported NOT BUILT.
3. Every figure also builds with all gates open (into a temporary folder).
4. plots.py and captions.py read only analysis outputs (figures/derived/), apart from the
   report HTML and the cached render; analysis.py reads nothing outside paper/.
5. No script holds an absolute or parent data path.
"""
from __future__ import annotations

import os
import re
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import yaml  # noqa: E402

import analysis  # noqa: E402
import common  # noqa: E402

SCRIPTS = ["analysis", "plots", "captions", "common", "build_all"]
GATED = {"5A", "5B", "5C", "5D", "S4A", "S4B", "S5A", "S5B"}
PLOT_INPUTS_ALLOWED = ("figures/derived/", "figures/rendered/", "09_report/report_pos.html")


def _build_figures(out: Path, gates_file: Path | None = None):
    os.environ["FIG_OUT_DIR"] = str(out)
    if gates_file:
        os.environ["FIG_GATES_FILE"] = str(gates_file)
    try:
        import plots
        common.OUT_DIR = out
        common.reset_reads()
        rows = []
        for build in plots.FIGURES:
            paths, r = build()
            assert all(p.exists() and p.stat().st_size > 0 for p in paths), build.__name__
            rows += r
        return rows, common.reads()
    finally:
        os.environ.pop("FIG_GATES_FILE", None)
        os.environ.pop("FIG_OUT_DIR", None)


def test_analysis_reproduces_expected_numbers():
    analysis.run()
    assert analysis.check() == [], analysis.check()
    paper = common.PAPER.resolve()
    for key, files in analysis.SOURCES.items():
        for f in files:
            assert paper in (paper / f).resolve().parents, f"{key} reads outside paper/: {f}"


def test_plots_read_only_analysis_outputs():
    analysis.run()
    with tempfile.TemporaryDirectory() as d:
        _, reads = _build_figures(Path(d))
    for panel, files in reads.items():
        for f in files:
            assert f.startswith(PLOT_INPUTS_ALLOWED), f"{panel} reads {f}, not an analysis output"


def test_captions_read_only_analysis_outputs():
    text = (HERE / "captions.py").read_text()
    assert "read_table(" not in text and "read_yaml(" not in text
    assert re.search(r"DERIVED\s*=\s*FIG_DIR / \"derived\"", text)


def test_build_with_committed_gates():
    analysis.run()
    with tempfile.TemporaryDirectory() as d:
        rows, _ = _build_figures(Path(d))
    g = common.gates()
    status = {f"{r['figure']}{r['panel']}": r["status"] for r in rows}
    if not g["reconstitution_volume_confirmed"]:
        assert status["S4A"] == status["S4B"] == "NOT BUILT"
    if not g["stage5_fix_confirmed"]:
        for k in ("5A", "5B", "5C", "S5A", "S5B"):
            assert status[k] == "NOT BUILT", k
    for k in ("1A", "1B", "3A", "3B", "3C", "3D", "4A", "4B", "4C"):
        assert status[k] == "REAL", k


def test_build_with_all_gates_open():
    analysis.run()
    lipid = str(common.read_table("figures/derived/app_trace.tsv").lipid_name.iloc[0])
    with tempfile.TemporaryDirectory() as d:
        gf = Path(d) / "gates.yaml"
        gf.write_text(yaml.safe_dump({"reconstitution_volume_confirmed": True,
                                      "stage5_fix_confirmed": True, "fig7E_lipid": lipid}))
        rows, _ = _build_figures(Path(d), gf)
    status = {f"{r['figure']}{r['panel']}": r["status"] for r in rows}
    for k in GATED:
        assert status[k] == "REAL", k


def test_no_absolute_or_parent_paths_in_scripts():
    bad = re.compile(r"""["'](/Users/|/lisc/|~/|\.\./)""")
    for name in SCRIPTS:
        for i, line in enumerate((HERE / f"{name}.py").read_text().splitlines(), 1):
            if "Chrome" in line:  # the renderer binary, not a data input
                continue
            assert not bad.search(line), f"{name}.py:{i}: {line.strip()}"


def test_read_outside_paper_is_refused():
    try:
        common.read_table("../lib/soillifeatlas/__init__.py")
    except (PermissionError, FileNotFoundError) as e:
        assert isinstance(e, PermissionError), e
    else:
        raise AssertionError("read outside paper/ was not refused")


if __name__ == "__main__":
    tests = [v for k, v in dict(globals()).items() if k.startswith("test_")]
    for t in tests:
        t()
        print("ok", t.__name__)
    print(f"{len(tests)} tests passed")
