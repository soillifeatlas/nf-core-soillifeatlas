"""Unit test for bin/fetch_atlas.py.

TDD: written before the wrapper. Exercises three of the four branches
end-to-end via subprocess:

  1. --atlas-path (offline/local override)
  2. Cache hit (pre-populated cache directory)
  3. Missing required artefact (error surface)
  4. --help (no-crash sanity)

The Zenodo download branch is *not* tested here because we do not yet
have a real atlas DOI; that branch is validated manually when the v0.1.0
Zenodo release is cut.
"""
import json
import subprocess
import sys
from pathlib import Path

import pandas as pd


REPO_ROOT = Path(__file__).resolve().parents[2]
BIN = REPO_ROOT / "bin" / "fetch_atlas.py"


def _build_tiny_atlas_dir(root: Path) -> Path:
    """Construct a minimal atlas directory that satisfies the required-artefact check."""
    atlas = root / "atlas_v0.1.0"
    atlas.mkdir(parents=True)
    pd.DataFrame(
        {
            "feature_id": ["f0", "f1"],
            "phylum": ["Bacillota", "Ascomycota"],
            "kingdom": ["Bacteria", "Fungi"],
        }
    ).to_parquet(atlas / "simper_fingerprint_atlas.parquet")
    return atlas


def test_fetch_atlas_local_override(tmp_path):
    """--atlas-path should bypass download and symlink/copy into --output-dir."""
    atlas = _build_tiny_atlas_dir(tmp_path / "source")
    output_dir = tmp_path / "out"
    cache = tmp_path / "cache"

    subprocess.run(
        [
            sys.executable,
            str(BIN),
            "--atlas-version",
            "v0.1.0",
            "--atlas-doi",
            "10.5281/zenodo.999999",
            "--cache-dir",
            str(cache),
            "--atlas-path",
            str(atlas),
            "--output-dir",
            str(output_dir),
        ],
        check=True,
    )

    # The required artefact must be reachable via the output dir.
    assert (output_dir / "simper_fingerprint_atlas.parquet").exists()

    manifest = json.loads((output_dir / "fetch_manifest.json").read_text())
    assert manifest["source"] == "local"
    assert manifest["atlas_version"] == "v0.1.0"
    assert "simper_fingerprint_atlas.parquet" in manifest["artefacts"]


def test_fetch_atlas_cache_hit(tmp_path):
    """A pre-populated cache dir at $CACHE/atlas_<version>/ should be used without download."""
    cache = tmp_path / "cache"
    # _build_tiny_atlas_dir creates <cache>/atlas_v0.1.0/ directly, which is the cache layout.
    _build_tiny_atlas_dir(cache)
    output_dir = tmp_path / "out"

    subprocess.run(
        [
            sys.executable,
            str(BIN),
            "--atlas-version",
            "v0.1.0",
            "--atlas-doi",
            "10.5281/zenodo.999999",
            "--cache-dir",
            str(cache),
            "--output-dir",
            str(output_dir),
        ],
        check=True,
    )

    assert (output_dir / "simper_fingerprint_atlas.parquet").exists()
    manifest = json.loads((output_dir / "fetch_manifest.json").read_text())
    assert manifest["source"] == "cache"
    assert manifest["atlas_version"] == "v0.1.0"


def test_fetch_atlas_missing_required_artefact(tmp_path):
    """A local override dir missing the required artefact must fail loudly."""
    empty = tmp_path / "empty_atlas"
    empty.mkdir()
    output_dir = tmp_path / "out"
    cache = tmp_path / "cache"

    result = subprocess.run(
        [
            sys.executable,
            str(BIN),
            "--atlas-version",
            "v0.1.0",
            "--atlas-doi",
            "10.5281/zenodo.999999",
            "--cache-dir",
            str(cache),
            "--atlas-path",
            str(empty),
            "--output-dir",
            str(output_dir),
        ],
        capture_output=True,
        text=True,
    )

    assert result.returncode != 0
    combined = (result.stderr + result.stdout).lower()
    assert "missing" in combined or "simper" in combined


def test_fetch_atlas_help_does_not_crash():
    """--help must exit 0 so Nextflow module-lint can introspect the script."""
    subprocess.run([sys.executable, str(BIN), "--help"], check=True)
