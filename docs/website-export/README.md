# Website handoff artefacts

This directory contains standalone artefacts for consumption by the
soillifeatlas.org website session (separate repo).

## Files

- `pipeline_blurb.md` — 200-word body copy for Stage 7 block
- `badges.md` — shields.io badge URLs (copy into Stage 7 component)
- `metrics.json` — JSON consumable by website JS for numeric callouts
- `dag.svg` — placeholder DAG figure (APPLY path; replace with final figure)
- `dag.README.md` — figure spec + swap instructions (source of truth is
  [`docs/ARCHITECTURE.md §1`](../ARCHITECTURE.md))

## Sync strategy (v0.1)

Manual copy from this directory into `website-repo/src/content/pipeline/`
at each pipeline release. Website rebuild picks up new content.

Future (v0.2+): Astro build-time fetch from raw.githubusercontent.com.
