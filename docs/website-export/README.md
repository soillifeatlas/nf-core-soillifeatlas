# Website handoff artefacts

Standalone artefacts consumed by the `pipeline.soillifeatlas.org`
subdomain. This directory is the single source of truth; the website
build step pulls these files at build time (raw.githubusercontent.com in
v0.2; manual copy in v0.1).

## Files

### Content for the Stage 7 block (pre-existing)

- `pipeline_blurb.md` — 200-word body copy
- `badges.md` — shields.io badge URLs
- `metrics.json` — numeric callouts (kingdom proportions, BC plausibility, feature counts)
- `dag.svg` — first-pass APPLY DAG (kept for backwards compatibility; superseded by `pipeline_dag_full.svg`)
- `dag.README.md` — spec for the original DAG figure

### New for `pipeline.soillifeatlas.org` (v0.1.0-dev artefacts)

- `form_schema.json` — JSON Schema draft 2020-12 defining the
  multi-step submission form (upload · instrument · advanced). Every
  property carries an `x-nextflow-param` hint that maps 1:1 to a param
  in `nextflow.config` / `nextflow_schema.json`. The website renders
  this schema into the form; validation uses the same schema.
- `pipeline_dag_simplified.svg` — 5-stage hero diagram
  (Ingest → Match → Correct → Decompose → Report). Horizontal,
  1200×240, scales cleanly to 600px. Headline visual for the landing
  hero section.
- `pipeline_dag_full.svg` — Per-process technical DAG for the
  "under the hood" section. Shows INGEST, PREPROC subworkflow,
  FETCH_ATLAS (moss side-input), the DECOMPOSE_APPLY subworkflow with
  all three correction layers, the 4 parallel decomposition methods,
  the three evaluation modules, and REPORT. ANNOTATE is dimmed as an
  optional parallel branch.
- `command_generator.js` — Pure, dependency-free function
  `generateNextflowCommand(formState) → string`. Works as an ES module
  (`import`) and via CommonJS (`require`). Only emits `--flag value`
  pairs for values that differ from pipeline defaults.
- `results_mock.html` — Standalone HTML results preview. No external
  dependencies; inline CSS with design tokens hard-coded. Template for
  the website's `/results-demo` route.
- `demo_bundle_manifest.md` — Spec for the
  `demo_climgrass_12samples_v0.1.0.tar.gz` bundle (ClimGrass 12-sample
  inputs + atlas artefacts + expected outputs). Tarball assembly is a
  separate task; this file is the build brief.

## Design tokens

All visual artefacts use the tokens declared in
`docs/plans/2026-04-17-soil-life-atlas-website-design.md`:

| Token | Hex |
|---|---|
| ink | `#1a1a1f` |
| paper | `#faf8f4` |
| soil | `#4a3728` |
| moss | `#2f4a3a` |
| spectral | `#2856d6` |
| muted | `#6b6860` |

Fonts: Fraunces (display), Inter (body), JetBrains Mono (technical /
parameter names). Files declare these `font-family`s but fall back to
generic system stacks so they render standalone.

## Sync strategy

v0.1: manual copy from this directory into the website repo at each
pipeline release. Website rebuild picks up new content.

v0.2+: Astro build-time fetch from raw.githubusercontent.com pinned to
the current pipeline tag.

## Parameter-name contract

`form_schema.json` and `command_generator.js` reference Nextflow param
names verbatim from `nextflow.config` (e.g. `IS_reference_compound`,
`RIE_floor`, `restrict_archaea_to_archlips`, `decomposition_methods`,
`primary_method`, `skip_gnps`, `skip_annotate`, `atlas_version`,
`soil_intensity`, `soil_mgf`, `sample_metadata`). If you rename a param
in the pipeline, update both files in lockstep or the website form will
silently generate broken commands.
