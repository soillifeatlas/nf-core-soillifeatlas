"""Write paper/<issue>/provenance.yaml next to an export."""
from __future__ import annotations

import hashlib
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path


def sha256(path: str | Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def git(*args: str) -> str:
    repo = Path(__file__).resolve().parents[2]
    return subprocess.run(["git", "-C", str(repo), *args], capture_output=True, text=True).stdout.strip()


def write(outdir: str | Path, inputs: dict[str, str | Path], params: dict, notes: str = "",
          pipeline_run: dict | None = None) -> Path:
    """Record SHA, dirty state, command, params and input checksums (YAML by hand, no dependency)."""
    out = Path(outdir) / "provenance.yaml"
    lines = [
        f"date: {datetime.now(timezone.utc).isoformat(timespec='seconds')}",
        f"export_git_sha: {git('rev-parse', 'HEAD')}",
        f"export_git_dirty: {bool(git('status', '--porcelain', '--', 'paper/scripts', 'lib'))}",
        f"command: {' '.join(sys.argv)!r}",
    ]
    if pipeline_run:
        lines.append("pipeline_run:")
        lines += [f"  {k}: {v!r}" for k, v in pipeline_run.items()]
    lines.append("params:")
    lines += [f"  {k}: {v!r}" for k, v in params.items()]
    lines.append("inputs:")
    for name, p in inputs.items():
        lines += [f"  {name}:", f"    path: {str(p)!r}", f"    sha256: {sha256(p)}"]
    if notes:
        lines.append(f"notes: {notes!r}")
    out.write_text("\n".join(lines) + "\n")
    return out


def read_run_provenance(path: str | Path) -> dict:
    """Pick commit/branch/start/finish from a pipeline run_provenance.txt."""
    keep = {}
    for line in Path(path).read_text().splitlines():
        for key in ("commit", "branch", "started", "finished", "uncommitted changes in clone", "nextflow"):
            if line.startswith(key + ":"):
                keep[key.replace(" ", "_")] = line.split(":", 1)[1].strip()
    keep["run_provenance_file"] = str(path)
    return keep
