# DAG figure — placeholder + swap instructions

The `dag.svg` in this directory is a hand-written placeholder rendering of
the APPLY workflow critical path (v0.1 live path). It is intentionally
minimal; the website session is expected to replace it with a polished
"tube map" style figure before the public launch.

## What the figure must show

1. **Input:** `samplesheet.csv` (raw_path | ion_mode | batch | role)
2. **PREPROC subworkflow:** INGEST → MZMINE_PP → MOLNET → ALIGN_LOESS,
   emitting aligned_feature_table + consensus MGF.
3. **FETCH_ATLAS** (side-input): Zenodo DOI + SHA256 verify, cached in
   `NXF_CACHE_DIR`.
4. **SIMPER_MATCH:** soil consensus × SIMPER atlas (matchms, cos ≥ 0.7,
   prec ≤ 5 ppm).
5. **DECOMPOSE_APPLY subworkflow:** L2 IS → L3 RIE(floor=0.20) → L5
   Archaea ref-filter → DECOMPOSE × 4 methods (nnls, std_bc,
   enriched_bc, fc_weighted_bc; primary = fc_weighted_bc) →
   EVAL_PLAUSIBILITY + EVAL_TREATMENT + DIAGNOSTIC_TOP_FEATURES.
6. **ANNOTATE subworkflow** (dashed, parallel, `--skip_annotate` default):
   MS-DIAL · SIRIUS · diag-ion · .msp match · ArchLips.
7. **REPORT:** MultiQC + custom panels, emitting
   `soillifeatlas_report.html` + `provenance.yaml`.

## Colour legend (suggested)

- **Blue** — v0.1 live APPLY modules
- **Amber** — atlas fetch / reference material
- **Green** — report / provenance
- **Grey dashed** — optional subworkflow (annotation)

## Source ASCII art

The canonical text representation is in
[`docs/ARCHITECTURE.md §1`](../ARCHITECTURE.md#1-pipeline-dag-simper-first-v01-apply-path-in-bold)
and [`docs/plans/2026-04-17-nextflow-pipeline-design.md §1`](../../../docs/plans/2026-04-17-nextflow-pipeline-design.md).
Treat those as source of truth; the SVG/PNG in this directory must stay
in sync at each release.

## PNG placeholder

`dag.png` is not committed in v0.1. The website session should rasterise
the final SVG to PNG at 2× density for retina rendering (`rsvg-convert -d 144
dag.svg -o dag.png` or equivalent).
