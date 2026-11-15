# soillifeatlas/nf-core-soillifeatlas: Usage

> Every parameter, with its default, is listed in `nextflow_schema.json`
> (`nextflow run . --help --show_hidden`).

## What the pipeline does

Apply mode profiles a soil LC–MS/MS lipid dataset against the Soil Life Atlas
reference. It reports, per soil sample, the share of lipid signal attributed to
each organism group (kingdom) and phylum, with four decomposition methods side by
side, plus plausibility and diagnostic tables.

## Which level of inference your data can reach

The atlas supports three separate levels of inference. Check yours before
running the pipeline.

| Level                 | Question                                            | Needs                                                                                                                                                   | Pipeline output                                             |
| --------------------- | --------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------- | ----------------------------------------------------------- |
| 1. Detection          | Is signal from a group present?                     | Any LC–MS/MS lipid data with MS/MS                                                                                                                      | `simper/*_verified_matches.parquet` (feature-level matches) |
| 2. Lipid-signal share | How much of the lipid signal comes from each group? | Data from the atlas workflow: modified monophasic Bligh–Dyer extraction with EquiSPLASH, Accucore C30, Q Exactive DDA, aligned with the atlas consensus | `decompose/composition_{kingdom,phylum}_<method>.parquet`   |
| 3. Relative biomass   | How much biomass does each group represent?         | Level 2 plus lipid-content conversion; assumes lipids come from living organisms                                                                        | Not produced by this version                                |

Decomposition (level 2) needs feature IDs in `--soil_intensity` that come from
the same alignment as the atlas consensus. Data from another extraction or
instrument should stop at level 1.

## Inputs

| Parameter           | Format    | Content                                                                |
| ------------------- | --------- | ---------------------------------------------------------------------- |
| `--soil_intensity`  | Parquet   | One row per `feature_id` (atlas-aligned), one column per sample        |
| `--soil_mgf`        | MGF       | MS/MS spectra for the same features                                    |
| `--sample_metadata` | TSV       | `sample_id` plus one column per contrast factor with values `Yes`/`No` |
| `--outdir`          | directory | Where results are written                                              |

The atlas bundle is fetched automatically (`--atlas_path` points at a local copy
for offline runs).

## Running

Check that everything works on the bundled tiny fixtures (a few minutes on a
laptop; host Python with `pip install -e .`):

```bash
nextflow run soillifeatlas/nf-core-soillifeatlas -profile test --outdir results
```

Your own data:

```bash
nextflow run soillifeatlas/nf-core-soillifeatlas \
  -profile local \
  --soil_intensity soil_quant.parquet \
  --soil_mgf soil.mgf \
  --sample_metadata sample_metadata.tsv \
  --outdir results
```

`-profile local` runs on the host Python environment. No module declares a
container yet, so `-profile docker` would also use the host environment; use
`local` until a container image is published.

## Normalisation

Intensities are divided by the internal standard named in
`--IS_reference_compound` (default `PC_15_0_18_1_d7`, the EquiSPLASH component
the paper normalises to). The relative-ionisation-efficiency correction (L3) is
experimental and off by default; `--enable_rie` turns it on.

## Treatment table

`eval/treatment_effects.tsv` compares groups named in `--treatment_contrast`. It
is a method-sensitivity check, not evidence for treatment effects.
