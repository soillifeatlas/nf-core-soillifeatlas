"""Checks for the figure scripts.  Run:  python paper/figures/test_figures.py  (or pytest).

1. Every figure builds with the gates as committed; gated figures are reported NOT BUILT.
2. Every figure also builds with all gates open (into a temporary folder), so the waiting
   figures are ready to re-run once the user confirms the inputs.
3. No panel reads a file outside paper/, and the scripts contain no absolute data paths.
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

import common  # noqa: E402

SCRIPTS = ["figures", "panels_atlas", "panels_workflow", "panels_qc", "panels_annotation",
           "panels_quant", "panels_interp", "panels_application"]
GATED = {"5A", "5B", "5C", "5D", "S4A", "S4B", "S5A", "S5B"}


def _build_all(out: Path, gates_file: Path | None = None):
    os.environ["FIG_OUT_DIR"] = str(out)
    if gates_file:
        os.environ["FIG_GATES_FILE"] = str(gates_file)
    try:
        import figures
        common.OUT_DIR = out
        common.reset_reads()
        rows = []
        for build in figures.FIGURES:
            paths, r = build()
            assert all(p.exists() and p.stat().st_size > 0 for p in paths), build.__name__
            rows += r
        return rows, common.reads()
    finally:
        os.environ.pop("FIG_GATES_FILE", None)
        os.environ.pop("FIG_OUT_DIR", None)


def _assert_inside_paper(reads):
    paper = common.PAPER.resolve()
    for panel, files in reads.items():
        for f in files:
            p = (paper / f).resolve()
            assert paper in p.parents, f"{panel} reads outside paper/: {f}"


def test_build_with_committed_gates():
    with tempfile.TemporaryDirectory() as d:
        rows, reads = _build_all(Path(d))
    _assert_inside_paper(reads)
    g = common.gates()
    status = {f"{r['figure']}{r['panel']}": r["status"] for r in rows}
    assert "PENDING" not in status.values()
    if not g["reconstitution_volume_confirmed"]:
        assert status["S4A"] == status["S4B"] == "NOT BUILT"
    if not g["stage5_fix_confirmed"]:
        for k in ("5A", "5B", "5C", "S5A", "S5B"):
            assert status[k] == "NOT BUILT", k
    for k in ("3A", "3B", "3C", "3D", "4A", "4B", "4C"):
        assert status[k] == "REAL", k


def test_build_with_all_gates_open():
    with tempfile.TemporaryDirectory() as d:
        gf = Path(d) / "gates.yaml"
        gf.write_text(yaml.safe_dump({"reconstitution_volume_confirmed": True,
                                      "stage5_fix_confirmed": True,
                                      "fig7E_lipid": _first_responsive_lipid()}))
        rows, reads = _build_all(Path(d), gf)
    _assert_inside_paper(reads)
    status = {f"{r['figure']}{r['panel']}": r["status"] for r in rows}
    for k in GATED:
        assert status[k] == "REAL", k


def _first_responsive_lipid() -> str:
    common.set_panel("")
    return str(common.read_table("08_application/pos/responsive.tsv").lipid_name.iloc[0])


def test_no_absolute_or_parent_paths_in_scripts():
    bad = re.compile(r"""["'](/Users/|/lisc/|~/|\.\./)""")
    for name in SCRIPTS + ["common", "build_all"]:
        text = (HERE / f"{name}.py").read_text()
        for i, line in enumerate(text.splitlines(), 1):
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
