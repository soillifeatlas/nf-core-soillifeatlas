# `demo_climgrass_12samples_v0.1.0.tar.gz` — Contents

**Status:** documented here; **tarball not yet assembled**. This file is the
spec the assembly task works from; the binary bundle lands on GitHub Releases
against the `v0.1.0-dev` tag once the Zenodo atlas DOI mints.

**Target install:** the user downloads the tarball, extracts it locally, and
runs the pipeline end-to-end without needing network access (atlas is
bundled). Matches the `--atlas_path <dir>` offline-override code path in
`bin/fetch_atlas.py`.

**Expected total size (gzipped):** ~40–80 MB (dominated by the soil MGF plus
the atlas consensus MGF).

---

## Inputs — what the user provides for a real run

| File | Size | Source | Notes |
|---|---|---|---|
| `soil_intensity.parquet` | ~1.5 MB (CSV → parquet) | `analysis-17/.../soil-data-without-background_quant.csv` | 12 samples × ~9k features, MZmine peak-area matrix. Converted CSV → parquet at bundle build. Expected columns: `feature_id` (str), one column per sample (float, peak area). Feature IDs must match the MGF scan IDs 1:1. |
| `soil.mgf` | ~8 MB | `analysis-17/.../soil-data-without-background.mgf` | GNPS-style MGF with per-feature consensus MS2 spectra. Scan/FEATURE_ID fields mirror the intensity matrix rows. |
| `sample_metadata.tsv` | ~1 KB | `analysis-17/.../sample_metadata.csv` (converted CSV → TSV) | 12 rows, header + 12 data rows. Columns: `sample_id, column_name, climate, drought, treatment, weight_g, notes`. The `treatment` column drives `EVAL_TREATMENT`. |

## Atlas artefacts — fetched by `FETCH_ATLAS` from Zenodo in online runs; bundled here for offline demo

These are the canonical files the `APPLY` workflow reads via the atlas
directory channel (see `workflows/apply.nf` step 3). Names must match the
bundle layout produced by `scripts/build_tiny_atlas_bundle.py`.

| File | Size | Purpose |
|---|---|---|
| `atlas/simper_fingerprint_atlas.parquet` | ~3–5 MB | 12,710 SIMPER fingerprint features × 18 phyla. The core reference matrix for decomposition. |
| `atlas/biomarker_atlas.parquet` | ~2 MB | 10,454 composite + IndVal biomarkers. Not consumed by decomposition in v0.1; shipped for external validation. |
| `atlas/rie_table_s10.csv` | ~15 KB | 109 class × adduct entries for the L3 RIE correction (Table S10 reference calibration). |
| `atlas/equisplash_IS_masses_POS.csv` | ~2 KB | 30 rows, one per (compound × adduct). 13 deuterated EquiSPLASH standards total (PC, LPC, PE, LPE, PG, PI, PS, TG, DG, MG, CE, SM, Cer). Defaults to the `LPE_18d7 [M+H]+` primary IS feature for v0.1. |
| `atlas/archlips_validated_features.csv` | ~20 KB | Feature-ID list for the L5 Archaea reference-quality filter. |
| `atlas/expected_kingdom_composition.csv` | ~1 KB | Reference kingdom proportions used by `EVAL_PLAUSIBILITY` to compute BC distance. |
| `atlas/consensus_unified_annotations.csv` | ~1–2 MB | Feature → class / adduct lookup consumed by `CORRECT_L3_RIE`. |
| `atlas/atlas_consensus.mgf` | ~10–20 MB | Consensus MS2 library for the SIMPER matcher. |
| `atlas/consensus_aligned_table.parquet` | ~3 MB | Atlas-wide aligned feature × sample intensity table. |
| `atlas/sample_phylum_map.csv` | ~5 KB | Atlas sample_id → phylum label mapping used by `DECOMPOSE`. |

## Expected outputs — what the user gets back after running

These are the artefacts emitted under `--outdir` after a successful
end-to-end run. Numbers below are the pinned ClimGrass reference values
asserted by the nightly full-demo CI job.

| File | Contents |
|---|---|
| `composition_fc_weighted_bc.parquet` | **Primary output.** 12 samples × 6 kingdoms × 18 phyla. Expected headline: Bacteria 35.7 ± 1.71, Fungi 29.6 ± 2.81, Plantae 19.0 ± 1.60, Archaea 9.1 ± 2.62, Animalia 3.0 ± 0.86, Protozoa 3.5 ± 0.49. |
| `composition_nnls.parquet` | Same shape, NNLS decomposition (sensitivity check). |
| `composition_std_bc.parquet` | Same shape, standard Bray-Curtis decomposition. |
| `composition_enriched_bc.parquet` | Same shape, IndVal-weighted Bray-Curtis. |
| `plausibility.tsv` | BC distance vs expected, one row per decomposition method. `fc_weighted_bc → 0.131`, below the 0.15 threshold. |
| `treatment_effects.tsv` | Mann-Whitney U per phylum × contrast. Expected direction: Actinomycetota ↑ drought (preserved from the published ClimGrass analysis); Basidiomycota ↓ drought (amplified vs raw); Bacteria ↓ warming, Fungi ↑ warming (newly revealed after L3 RIE correction). |
| `diagnostic_top_features.tsv` | Top-10 features per phylum (RIE over-amp sanity check; no single feature >3.1% of any phylum). |
| `multiqc_report.html` | Rendered MultiQC + custom panels. |
| `pipeline_info/provenance.yaml` | Git SHA, atlas DOI + SHA256, container digests, params, input hashes, runtime. |
| `pipeline_info/execution_report_*.html` | Standard Nextflow trace + DAG + timeline. |

---

## Assembly plan

1. **Copy soil inputs** from the ClimGrass working directory:
   ```
   SRC=/Users/rahulsamrat/Desktop/Projects/soilmass-analysis/analysis/analysis-17/positive/soil_decomposition/climgrass-experiment
   cp $SRC/soil-data-without-background_quant.csv   staging/inputs/
   cp $SRC/soil-data-without-background.mgf         staging/inputs/soil.mgf
   cp $SRC/sample_metadata.csv                      staging/inputs/
   ```
   Then convert CSV → parquet and CSV → TSV using the existing
   `scripts/*_to_parquet.py` helpers (or polars in an ad-hoc script).

   > **Path note:** the task spec points at `soilmass-viewer/...`; the
   > actual data lives under `soilmass-analysis/...` on this machine.
   > The layout underneath is identical.

2. **Copy atlas artefacts** from the build locations under
   `soilmass-analysis/analysis/analysis-15/` (SIMPER fingerprint, biomarker,
   consensus MGF, consensus aligned table, sample_phylum_map) and
   `soilmass-analysis/analysis/analysis-19/00_inputs/` (rie_table_s10.csv,
   equisplash_IS_masses_POS.csv, archlips_validated_features.csv,
   expected_kingdom_composition.csv, consensus_unified_annotations.csv).
   Stage them under `staging/atlas/` with the canonical names above.

3. **Dry-run the pipeline** against the staged bundle:
   ```
   nextflow run soillifeatlas/nf-core-soillifeatlas \
     -r v0.1.0-dev -profile local \
     --mode apply \
     --atlas_path staging/atlas \
     --soil_intensity staging/inputs/soil_intensity.parquet \
     --soil_mgf staging/inputs/soil.mgf \
     --sample_metadata staging/inputs/sample_metadata.tsv \
     --outdir staging/results
   ```
   Assert: all 4 compositions emitted, `plausibility.tsv` has
   `fc_weighted_bc → 0.131 ± 0.005`, `treatment_effects.tsv` shows
   Actinomycetota ↑ drought direction.

4. **Tarball + attach to GitHub Releases**:
   ```
   tar -czf demo_climgrass_12samples_v0.1.0.tar.gz -C staging \
     inputs atlas README.md
   gh release upload v0.1.0-dev demo_climgrass_12samples_v0.1.0.tar.gz
   ```

5. **Verify the tarball** by running the pipeline from a clean directory:
   extract, `nextflow run ... --atlas_path ./atlas --soil_intensity ./inputs/soil_intensity.parquet ...`
   and re-compare the output composition against the pinned numbers.

---

## Source file cross-check

Confirmed on this machine:

- `/Users/rahulsamrat/Desktop/Projects/soilmass-analysis/analysis/analysis-17/positive/soil_decomposition/climgrass-experiment/soil-data-without-background_quant.csv` (1.56 MB)
- `/Users/rahulsamrat/Desktop/Projects/soilmass-analysis/analysis/analysis-17/positive/soil_decomposition/climgrass-experiment/soil-data-without-background.mgf` (8.18 MB)
- `/Users/rahulsamrat/Desktop/Projects/soilmass-analysis/analysis/analysis-17/positive/soil_decomposition/climgrass-experiment/sample_metadata.csv` (897 B, 13 lines = header + 12 samples)
- `/Users/rahulsamrat/Desktop/Projects/soilmass-analysis/analysis/analysis-19/00_inputs/equisplash_IS_masses_POS.csv` (31 lines, 13 compounds × 2–3 adducts each)

Atlas artefact sizes are unverified and reported as ranges until real files
are staged; replace the ranges with measured values once the staging copy
lands.
