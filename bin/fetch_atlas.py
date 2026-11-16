#!/usr/bin/env python3
"""Fetch the soillifeatlas reference atlas with SHA256 verification.

Three resolution branches, in priority order:

  1. ``--atlas-path <dir>``: offline override. Verify the directory contains the
     required artefacts, then expose it via ``--output-dir``. No network, no cache.
     Used for air-gapped runs and CI/demo fixtures.
  2. Cache hit: if ``$CACHE_DIR/atlas_<VERSION>/`` already has the required
     artefacts, skip download and expose the cached copy. Subsequent pipeline
     runs on the same machine hit this branch.
  3. Zenodo download: resolve the DOI to a numeric record ID, fetch
     ``soillifeatlas_atlas_<VERSION>.tar.gz``, verify SHA256 if provided,
     extract into the cache, then expose via ``--output-dir``.

A JSON manifest (``fetch_manifest.json``) is written inside ``--output-dir`` and
echoed to stderr so Nextflow logs capture provenance (source, version, sha, path).

Standard library only: no new pipeline dependencies.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import sys
import tarfile
import urllib.request
from pathlib import Path


# Canonical atlas artefacts the downstream decomposition + SIMPER workflow needs.
# Keep this list conservative: adding an entry here will break all existing
# atlases until they are republished with the new file. For v0.1 we only
# require the SIMPER fingerprint; other reference tables are looked up
# opportunistically by consumer scripts.
REQUIRED_ARTEFACTS = [
    "simper_fingerprint_atlas.parquet",
]


def _verify_artefacts(atlas_dir: Path) -> list[str]:
    """Raise FileNotFoundError if any required artefact is missing.

    Returns the list of artefacts actually present (== REQUIRED_ARTEFACTS on success).
    """
    missing = [a for a in REQUIRED_ARTEFACTS if not (atlas_dir / a).exists()]
    if missing:
        raise FileNotFoundError(
            f"Atlas at {atlas_dir} missing required artefact(s): {missing}"
        )
    return [a for a in REQUIRED_ARTEFACTS if (atlas_dir / a).exists()]


def _sha256(path: Path) -> str:
    """Stream-hash a file in 1 MiB chunks to keep memory flat."""
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _emit_manifest(manifest: dict, output_dir: Path) -> None:
    """Write the manifest inside output_dir *and* echo it to stderr for Nextflow logs."""
    (output_dir / "fetch_manifest.json").write_text(json.dumps(manifest, indent=2))
    print(json.dumps(manifest), file=sys.stderr)


def _link_or_copy(src: Path, dst: Path) -> None:
    """Expose ``src`` at ``dst``. Prefer a symlink; fall back to copytree on OSError
    (e.g. Windows without symlink privilege, or cross-device links).

    If ``dst`` already exists (file, dir, or stale symlink), it is removed first.
    The ``src -> dst`` direction: ``dst`` is what the pipeline reads;
    ``src`` is the authoritative source (cache or user-provided path).
    """
    if dst.is_symlink() or dst.exists():
        if dst.is_symlink() or dst.is_file():
            dst.unlink()
        else:
            shutil.rmtree(dst)
    try:
        dst.symlink_to(src, target_is_directory=True)
    except OSError:
        shutil.copytree(src, dst)


def _parse_zenodo_record_id(doi: str) -> str:
    """Pull the numeric record id out of a Zenodo DOI string.

    Accepts the canonical ``10.5281/zenodo.XXXXXXX`` form; returns ``XXXXXXX``.
    """
    return doi.rsplit(".", 1)[-1]


def main() -> int:
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--atlas-version", required=True, help="e.g. v0.1.0")
    p.add_argument(
        "--atlas-doi",
        required=True,
        help="Zenodo DOI, e.g. 10.5281/zenodo.XXXXXXX",
    )
    p.add_argument(
        "--atlas-sha256",
        default="",
        help="Expected SHA256 of the .tar.gz (empty string skips the check).",
    )
    p.add_argument(
        "--cache-dir",
        required=True,
        help="Persistent cache root; atlas lands at $CACHE_DIR/atlas_<VERSION>/.",
    )
    p.add_argument(
        "--atlas-path",
        default=None,
        help="Local dir to use instead of downloading (air-gapped / CI override).",
    )
    p.add_argument(
        "--output-dir",
        required=True,
        help="Where downstream Nextflow processes read the atlas from.",
    )
    args = p.parse_args()

    output_dir = Path(args.output_dir)
    output_dir.parent.mkdir(parents=True, exist_ok=True)

    # Branch 1: offline local override.
    if args.atlas_path:
        atlas_path = Path(args.atlas_path).resolve()
        if not atlas_path.is_dir():
            print(
                f"--atlas-path {atlas_path} is not a directory",
                file=sys.stderr,
            )
            return 1
        artefacts = _verify_artefacts(atlas_path)
        _link_or_copy(atlas_path, output_dir)
        _emit_manifest(
            {
                "source": "local",
                "atlas_version": args.atlas_version,
                "path": str(atlas_path),
                "artefacts": artefacts,
            },
            output_dir,
        )
        return 0

    cache_dir = Path(args.cache_dir).resolve()
    cache_dir.mkdir(parents=True, exist_ok=True)
    cached = cache_dir / f"atlas_{args.atlas_version}"

    # Branch 2: cache hit (cheap check — just look for the canonical artefact).
    if (cached / REQUIRED_ARTEFACTS[0]).exists():
        artefacts = _verify_artefacts(cached)
        _link_or_copy(cached, output_dir)
        _emit_manifest(
            {
                "source": "cache",
                "atlas_version": args.atlas_version,
                "path": str(cached),
                "artefacts": artefacts,
            },
            output_dir,
        )
        return 0

    # Branch 3: Zenodo download.
    record_id = _parse_zenodo_record_id(args.atlas_doi)
    url = (
        f"https://zenodo.org/record/{record_id}/files/"
        f"soillifeatlas_atlas_{args.atlas_version}.tar.gz"
    )
    tarball = cache_dir / f"atlas_{args.atlas_version}.tar.gz"
    print(f"fetch_atlas: downloading {url} -> {tarball}", file=sys.stderr)
    urllib.request.urlretrieve(url, tarball)

    observed_sha = _sha256(tarball)
    if args.atlas_sha256 and observed_sha != args.atlas_sha256:
        tarball.unlink()
        raise ValueError(
            f"SHA256 mismatch for {url}: expected {args.atlas_sha256}, "
            f"got {observed_sha}"
        )

    with tarfile.open(tarball) as t:
        t.extractall(cache_dir)

    artefacts = _verify_artefacts(cached)
    _link_or_copy(cached, output_dir)
    _emit_manifest(
        {
            "source": "zenodo",
            "atlas_version": args.atlas_version,
            "url": url,
            "sha256": observed_sha,
            "path": str(cached),
            "artefacts": artefacts,
        },
        output_dir,
    )
    return 0


if __name__ == "__main__":
    sys.exit(main() or 0)
