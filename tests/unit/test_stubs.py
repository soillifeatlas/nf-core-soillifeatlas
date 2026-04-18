"""Smoke tests for the v0.1 train-mode stubs.

``bin/align_loess.py``, ``bin/build_biomarker_atlas.py``, and
``bin/build_simper_fp_atlas.py`` are intentionally stubbed in v0.1: the
real implementations are deferred to v0.2. The Nextflow train-mode
modules still reference these scripts, so they must:

  1. Exist and be discoverable via ``--help`` (exit 0, clean argparse).
  2. Exit non-zero (distinct code 2) with an informative stderr message
     when actually invoked, so Nextflow surfaces a process failure rather
     than silently producing empty outputs.
"""
import subprocess
import sys
from pathlib import Path

import pytest


REPO_ROOT = Path(__file__).resolve().parents[2]

STUB_SCRIPTS = [
    "align_loess.py",
    "build_biomarker_atlas.py",
    "build_simper_fp_atlas.py",
]


@pytest.mark.parametrize("script_name", STUB_SCRIPTS)
def test_stub_help_works(script_name):
    """--help must exit 0 so Nextflow module-lint can introspect the script."""
    script = REPO_ROOT / "bin" / script_name
    subprocess.run([sys.executable, str(script), "--help"], check=True)


@pytest.mark.parametrize("script_name", STUB_SCRIPTS)
def test_stub_exits_non_zero_with_informative_message(script_name):
    """Running a stub with its required args still exits 2 with an explanation."""
    script = REPO_ROOT / "bin" / script_name

    # Pass enough args to satisfy argparse so we exercise the stub body, not argparse's
    # "missing required argument" exit. Values are placeholders - the stub body should
    # reject before anything looks at them.
    if script_name == "align_loess.py":
        extra = ["--inputs", "x", "--output", "/tmp/x"]
    elif script_name == "build_biomarker_atlas.py":
        extra = [
            "--consensus-aligned-table",
            "x",
            "--sample-metadata",
            "y",
            "--output-biomarker-atlas",
            "/tmp/x",
        ]
    else:
        extra = [
            "--consensus-aligned-table",
            "x",
            "--sample-metadata",
            "y",
            "--output-simper-fp-atlas",
            "/tmp/x",
        ]

    result = subprocess.run(
        [sys.executable, str(script), *extra],
        capture_output=True,
        text=True,
    )

    assert result.returncode == 2, (
        f"{script_name} should exit 2 for 'not implemented', got {result.returncode}. "
        f"stderr={result.stderr}"
    )
    msg = result.stderr.lower()
    assert (
        "v0.1 stub" in msg or "deferred" in msg or "not implemented" in msg
    ), f"{script_name} stderr did not explain stub status: {result.stderr}"
