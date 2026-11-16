# soillifeatlas/nf-core-soillifeatlas: Usage

> Every parameter, with its default, is listed in `nextflow_schema.json`
> (`nextflow run . --help --show_hidden`).

## What the pipeline does

The pipeline has two modes:

- **Process mode** (`--mode process`) takes raw Orbitrap files to quality-checked
  feature tables: conversion (stage 1), MZmine feature finding (stage 2) and a
  quality-check report (stage 3). Annotation, quantification with recovery, the
  report and the Atlas step (stages 4–7) are being added.
- **Apply mode** (`--mode apply`, the default) compares soil data with the Soil
  Life Atlas reference (see below).

## Process mode

### The sample sheet

One row per injection, both polarities in one sheet if you like:

```csv
file,polarity,role,group,dry_mass_g,extract_mg,spike_ng,volume_ul
/data/run1/S07_POS.raw,pos,sample,drought,2.01,4.3,10,200
/data/run1/QC_01_POS.raw,pos,pooled_qc,,,,,
/data/run1/NIST_01_POS.raw,pos,nist,,,,,
/data/run1/CAL_1_POS.raw,pos,standard_cal,,,,,
/data/run1/EBLK_01_POS.raw,pos,extraction_blank,,,,,
/data/run1/SBLK_01_POS.raw,pos,blank,,,,,
```

| Column                                              | Required | Meaning                                                                                                                                                     |
| --------------------------------------------------- | -------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `file`                                              | yes      | Thermo `.raw` (converted in stage 1) or `.mzML` (used as it is)                                                                                             |
| `polarity`                                          | yes      | `pos` or `neg`; each polarity is processed separately                                                                                                       |
| `role`                                              | yes      | `sample`, `pooled_qc`, `blank` (solvent), `extraction_blank`, `nist` (NIST SRM 1950), `standard_cal` (internal standard injected on its own, not extracted) |
| `group`, `batch`                                    | no       | treatment group and acquisition batch                                                                                                                       |
| `dry_mass_g`, `extract_mg`, `spike_ng`, `volume_ul` | no       | sample dry mass, dried lipid-fraction mass, internal standard spiked, redissolve volume; used for recovery and amounts per gram (stage 5)                   |

Injection order is read from the mzML run start times. Convert all files of a
run the same way: different converters can record the start time with different
time-zone offsets.

### Internal standards

`--is_targets` points to a table of the internal standards in your samples. The
default (`assets/equisplash_targets.csv`) lists all 13 EquiSPLASH components
with the common adducts in both polarities. Use `--is_compounds` to keep only
the ones you spiked (names from the `compound` column, comma-separated), or give
your own CSV with at least `compound`, `polarity`, `adduct` and `mz`. Optional
columns: `rt` (expected retention time on your method), `priority`
(`primary`/`secondary`), `lipid_class`, `spike_ng`.

Soil contains natural lipids within a few ppm of the deuterated standards, so a
standard is only accepted as confirmed when its retention time is known: from
an `rt` column (window `--is_rt_window`, default 0.3 min), or from
`standard_cal` injections in which the standard dominates. Otherwise the report
marks it as unconfirmed. `assets/examples/` has an example table with retention
times.

### Quality checks (stage 3)

Pass/fail follows Samrat et al. (2025): internal-standard retention time within
±0.1 min (`--qc_rt_tolerance`) and mass error within ±5 ppm
(`--qc_mass_accuracy_ppm`) in pooled QC, NIST and standard injections;
internal-standard peak-area RSD below 20 % in pooled QCs (`--qc_is_rsd`). The
report also gives the pooled-QC RSD of every feature, the share of sample
features also present in blanks, and features and total signal per injection in
run order. Variation of the standards across soils is reported for information:
it reflects matrix effects and recovery, which stage 5 corrects.

### MZmine

MZmine 4 needs a user login. Log in once on the machine that will run the
pipeline (`mzmine --login-console`, or through the MZmine window) and pass the
created file with `--mzmine_user ~/.mzmine/users/<name>.mzuser`. The file is
copied into the task folder for the run and not written to any output.
`--mzmine_batch_pos` and `--mzmine_batch_neg` replace the shipped settings.

User files expire after a while; MZmine then reports that the user "needs to
login again" and the pipeline stops with a message saying so. Log in again and
rerun with `-resume`: finished conversions are reused. On clusters with a small
home quota, write the user file to a work area instead:

```bash
JAVA_TOOL_OPTIONS="-Duser.home=/path/to/work/mzmine_home" mzmine --login-console
# user file: /path/to/work/mzmine_home/.mzmine/users/<name>.mzuser
```

### Running

```bash
nextflow run soillifeatlas/nf-core-soillifeatlas \
  --mode process \
  --input samplesheet.csv \
  --mzmine_user ~/.mzmine/users/me.mzuser \
  --run_name myrun \
  --outdir results
```

On LiSC use `-profile lisc` and a small personal config for site paths (see
`conf/lisc.config`).

Outputs: `results/mzmine/<run>_<pol>.features_full.csv` (and the GNPS FBMN, SIRIUS and
annotation exports), `results/qc/<run>_<pol>.qc_report.html` with the tables
behind it, and `results/mzml/` if `--save_mzml` is set.

## Apply mode

Apply mode profiles a soil LC–MS/MS lipid dataset against the Soil Life Atlas
reference. It reports, per soil sample, the share of lipid signal attributed to
each organism group (kingdom) and phylum, with four decomposition methods side by
side, plus plausibility and diagnostic tables.

### Which level of inference your data can reach

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

### Inputs

| Parameter           | Format    | Content                                                                |
| ------------------- | --------- | ---------------------------------------------------------------------- |
| `--soil_intensity`  | Parquet   | One row per `feature_id` (atlas-aligned), one column per sample        |
| `--soil_mgf`        | MGF       | MS/MS spectra for the same features                                    |
| `--sample_metadata` | TSV       | `sample_id` plus one column per contrast factor with values `Yes`/`No` |
| `--outdir`          | directory | Where results are written                                              |

The atlas bundle is fetched automatically (`--atlas_path` points at a local copy
for offline runs).

### Running

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

### Normalisation

Intensities are divided by the internal standard named in
`--IS_reference_compound` (default `PC_15_0_18_1_d7`, the EquiSPLASH component
the paper normalises to). The relative-ionisation-efficiency correction (L3) is
experimental and off by default; `--enable_rie` turns it on.

### Treatment table

`eval/treatment_effects.tsv` compares groups named in `--treatment_contrast`. It
is a method-sensitivity check, not evidence for treatment effects.
