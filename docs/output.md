# soillifeatlas/nf-core-soillifeatlas: Output

All paths are relative to `--outdir`. `<method>` is one of `nnls`, `std_bc`,
`enriched_bc`, `fc_weighted_bc` (primary: `--primary_method`, default
`fc_weighted_bc`).

## `fetch/`

- `atlas_fetched/` — the atlas bundle used for the run (consensus MGF, consensus
  intensity table, SIMPER fingerprint atlas, internal-standard table, annotation
  table, sample-to-phylum map, expected kingdom ranges).

## `simper/` — level 1, detection

- `demo_verified_matches.parquet` — soil features matched to atlas fingerprint
  features by MS/MS cosine (thresholds: `--simper_min_cos`,
  `--simper_min_matched_peaks`, `--simper_precursor_ppm`, `--simper_fragment_tol`).

## `correct/`

- `intensity_l2_is_corrected.parquet` — soil intensities divided by the internal
  standard.
- `intensity_l3_rie_corrected.parquet` — only with `--enable_rie` (experimental).
- `simper_atlas_l5_filtered.parquet` — fingerprint atlas after the reference-quality
  filter (Archaea restricted to ArchLips-validated features).

## `decompose/` — level 2, lipid-signal share

- `composition_kingdom_<method>.parquet` — long format: `sample_id`, `kingdom`,
  `proportion_pct`; proportions sum to 100 per sample.
- `composition_phylum_<method>.parquet` — the same at phylum level.

These are shares of lipid signal, not biomass.

## `eval/`

- `plausibility_<method>.tsv` — mean kingdom composition compared with literature
  ranges for grassland soil (`bc_vs_expected`, `in_range_fraction`, deviations).
- `treatment_effects.tsv` — Mann–Whitney U per phylum and contrast on the primary
  method. A method-sensitivity check, not a claim about treatment effects.

## `diagnostic/`

- `top_features_per_phylum.tsv` — the features contributing most signal to each
  phylum (`--diagnostic_top_n`, default 10).

## `pipeline_info/`

Nextflow execution report, timeline, trace and DAG, plus `params.json` and
software versions.
