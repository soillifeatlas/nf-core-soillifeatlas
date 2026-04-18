# soillifeatlas/nf-core-soillifeatlas — Architecture

**Status:** v0.1-dev. End-to-end APPLY workflow runs on tiny fixtures;
v0.1.0 release blocked on Zenodo DOI for the reference atlas.

This 1-page artefact is reviewer-facing: it maps the published method
(Samrat et al. 2025, *Soil Biol Biochem*, doi:10.1016/j.soilbio.2025.109892)
to the Nextflow modules that execute it, and documents the "not vaporware"
CI contract. Deeper design lives in the
[design doc](../../docs/plans/2026-04-17-nextflow-pipeline-design.md) and
[implementation plan](../../docs/plans/2026-04-17-nextflow-pipeline-implementation.md).

---

## 1. Pipeline DAG (SIMPER-first, v0.1 APPLY path in bold)

Decomposition depends **only on feature_id matching to the SIMPER fingerprint
atlas**, not on lipid annotations (no single feature contributes >3.1% of
any phylum's distinctiveness; signal is distributed + redundant). Consequence:
`ANNOTATE` is parallel enrichment, not a gate. Demo defaults to `--skip_annotate`
and runs in ~20 min on a laptop.

```
                  samplesheet.csv (raw_path | ion_mode | batch | role)
                              │
                              ▼
                    ┌──────────────────┐
                    │ PREPROC subwf    │   (shared by train + apply)
                    │  INGEST          │
                    │  MZMINE_PP       │
                    │  MOLNET          │
                    │  ALIGN_LOESS     │
                    └────────┬─────────┘
                             │ aligned_feature_table + consensus_MGF
        ┌────────────────────┼────────────────────┐
        ▼ TRAIN              │                    ▼ APPLY  *** v0.1 ***
 ┌──────┴──────┐       ┌─────┴──────┐      ┌─────┴────────────┐
 │ BIOMARKER   │       │ SIMPER_FP_ │      │ SIMPER_MATCH     │
 │  composite  │       │   ATLAS    │      │ (soil × atlas)   │
 │  IndVal.g   │       │ (12,710 ft │      └──────┬───────────┘
 │  tiers      │       │  × 18 ph)  │             │
 └──────┬──────┘       └─────┬──────┘             ▼
        │                    │            ┌───────────────────┐
        ▼                    │            │ DECOMPOSE_APPLY   │
 ┌─────────────┐             │            │  L2 IS            │
 │ VALIDATE_EXT│             │            │  L3 RIE(floor=.2) │
 │  fastMASST  │             │            │  L5 ref-filter    │
 │  DreaMS     │             │            │  4×DECOMPOSE      │
 │  Pan-ReDU   │             │            │  eval plaus+TE    │
 └──────┬──────┘             │            │  diagnostic       │
        │                    │            └────────┬──────────┘
        └─────────┬──────────┘                     │
                  ▼                                │
          ┌──────────────┐                         │
          │ ATLAS_RELEASE│ (Zenodo)                │
          │  simper_fp   │                         │
          │  biomarkers  │                         │
          │  msp libs    │                         │
          │  RIE / IS    │                         │
          └──────────────┘                         │
                                                   │
        ┌────── ANNOTATE ───────┐                  │
        │  MS-DIAL, SIRIUS,     │ (parallel,       │
        │  diag-ion, .msp,      │  non-blocking,   │
        │  ArchLips             │  `--skip` default)
        └──────────┬────────────┘                  │
                   │                               │
                   └──────────────┬────────────────┘
                                  ▼
                          ┌──────────────┐
                          │ REPORT       │
                          │  (MultiQC +  │
                          │   custom)    │
                          └──────────────┘
```

v0.1 executes the APPLY path end-to-end on tiny fixtures (right side). TRAIN
is scaffolded; real execution is an ~8 hr LISC job and ships in v0.2+.

---

## 2. SOP traceability — Samrat 2025 Methods → Nextflow modules → framework/*.py

Every step in the paper's Methods section maps to one Nextflow module, which
in turn wraps one function in `framework/*.py` (the reference Python
implementation, vendored to `bin/` at each release). The pipeline is a
reproducible container around the same code that produced the published
numbers.

| Samrat 2025 §Methods | `framework/` fn | Nextflow module | v0.1 status |
|---|---|---|---|
| Sample prep + LC-MS acquisition | — (wet lab) | — | out of scope |
| Feature detection (MZmine 4.9.14 IIMN) | — | `MZMINE_PP` | scaffolded |
| Molecular networking (GNPS2 FBMN + Classical MN) | — | `MOLNET` | scaffolded |
| Cross-batch alignment (3-stage LOESS) | `align_loess.py` | `ALIGN_LOESS` | scaffolded |
| Biomarker discovery (composite scoring) | `build_biomarker_atlas.py` | `BIOMARKER_COMPOSITE` | TRAIN v0.2+ |
| Biomarker discovery (IndVal.g cross-batch) | `build_biomarker_atlas.py` | `BIOMARKER_INDVAL` | TRAIN v0.2+ |
| SIMPER fingerprint atlas | `build_simper_fp_atlas.py` | `SIMPER_FP_ATLAS` | TRAIN v0.2+ |
| Annotation pipeline (SIRIUS, diag-ion, ArchLips, .msp, MS-DIAL) | various | `ANNOTATE` subworkflow | parallel, `--skip` |
| External validation (fastMASST, DreaMS, Pan-ReDU) | `validate_*.py` | `VALIDATE_EXT` subworkflow | TRAIN v0.2+ |
| Spectral matching (soil × SIMPER atlas) | `simper_match.py` | `SIMPER_MATCH` | **v0.1 live** |
| IS normalization (analysis-19 L2) | `corrections.apply_IS` | `CORRECT_L2_IS` | **v0.1 live** |
| RIE correction (analysis-19 L3, floor=0.20) | `corrections.apply_RIE` | `CORRECT_L3_RIE` | **v0.1 live** |
| Archaea reference-filter (analysis-19 L5) | `corrections.restrict_archaea_to_archlips` | `CORRECT_L5_REF_FILTER` | **v0.1 live** |
| 4 decomposition methods (nnls, std_bc, enriched_bc, fc_weighted_bc) | `decomposition.*` | `DECOMPOSE` (parallel via `each`) | **v0.1 live** |
| Plausibility scoring (BC vs expected) | `evaluation.plausibility` | `EVAL_PLAUSIBILITY` | **v0.1 live** |
| Treatment effect verification (Mann-Whitney U) | `evaluation.treatment_effects` | `EVAL_TREATMENT` | **v0.1 live** |
| Feature-domination diagnostic (RIE over-amp guard) | `evaluation.diagnose_top_features` | `DIAGNOSTIC_TOP_FEATURES` | **v0.1 live** |
| Atlas fetch (Zenodo DOI + SHA256 verify + cache) | `fetch_atlas.py` | `FETCH_ATLAS` | **v0.1 live** |

All reference constants (RIE floor=0.20, IS reference = LPE_18d7, atlas SHA,
IIMN XML digest, container digests) are pinned in `nextflow.config`. Override
requires explicit `--allow_param_override`.

---

## 3. Deployment matrix

| Profile | Executor | Container engine | Audience | Status |
|---|---|---|---|---|
| `local` | local | none (system Python) | Developer, CI, reviewer laptop | v0.1 live |
| `docker` | local | Docker | Reviewer laptop demo | v0.1 live |
| `lisc_slurm` | SLURM | Apptainer | Real compute (LISC at Uni Vienna; partitions `xeon_0384` default, `zen2_2048` heavy alignment) | v0.1 live |
| `google_batch` | Google Batch | Docker | Cloud-scale (post-funding) | v0.1 stub |

v0.1 CI executes the `local` profile on GitHub Actions `ubuntu-latest`. The
`google_batch` profile ships as a stub: config + docs exist, but execution
unlocks with Google.org funding (costs $0 today).

---

## 4. "Not vaporware" contract

Reviewer promise: green badges = pipeline reproduces Samrat 2025 *today*,
not at some future milestone.

**PR CI (`ci.yml`, every push + pull request):**

- Python unit tests: 32/32 green (correction layers, decomposition methods,
  evaluation metrics, atlas fetch)
- nf-test modules + subworkflows + workflow: 11/11 green
- APPLY workflow end-to-end on tiny fixtures in <5 min
- nf-core pipelines lint (known issues; `continue-on-error: true` in v0.1)

**Nightly full-demo (`full-demo.yml`, 03:00 UTC + on release tags):**

- APPLY workflow end-to-end on tiny fixtures, via `nf-test`
- Post-v0.1.0: 12-sample ClimGrass reproduction with kingdom-composition
  assertions (Bacteria 35.7% ±2pp; Fungi 29.6% ±3pp; Plantae 19.0% ±2pp;
  Archaea 9.1% ±3pp; Animalia 3.0% ±1pp; Protozoa 3.5% ±1pp; BC vs expected
  ≤0.15; Actinomycetota drought direction preserved; atlas SHA matches
  pinned digest; all 4 decomposition methods converged). Asserted via
  `tests/test_samrat_2025_reproduction.py` — scaffolded, enables on DOI.

**Provenance (`pipeline_info/provenance.yaml`, emitted every run):**

- git_sha, nextflow_version, pipeline_version
- atlas_version + DOI + SHA256
- container digests (all pinned `sha256:...`)
- params (RIE_floor, IS_reference_compound, IS_spiked_pmol, …)
- inputs (samplesheet_sha256, n_samples, ion_modes)
- run metadata (started, completed, executor, runtime_seconds)

Reviewer diff: reference `provenance.yaml` shipped in Zenodo atlas release
vs local run → identical compute verified.

---

## 5. Further reading

- **Design doc:**
  `/Users/rahulsamrat/Desktop/Projects/soilfoodwebatlas/docs/plans/2026-04-17-nextflow-pipeline-design.md`
  (full architecture, scope decisions, risks)
- **Implementation plan:**
  `/Users/rahulsamrat/Desktop/Projects/soilfoodwebatlas/docs/plans/2026-04-17-nextflow-pipeline-implementation.md`
  (phase-by-phase build log)
- **Paper:** Samrat et al. 2025, *Soil Biol Biochem*,
  [doi:10.1016/j.soilbio.2025.109892](https://doi.org/10.1016/j.soilbio.2025.109892)
- **Quantification correction feedback loop (analysis-19):** origin of L2 IS
  / L3 RIE(floor=0.20) / L5 ref-filter. See
  `/soilmass-viewer/analysis/analysis-19/REPORT.md`.
- **nf-core framework:** Ewels et al. 2020, *Nat Biotechnol*,
  [doi:10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x)
