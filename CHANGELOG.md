# soillifeatlas/soillifeatlas: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v0.1.0dev - unreleased

Initial development version, created with the [nf-core](https://nf-co.re/) template.

### `Added`

- `-profile test`: apply mode on the bundled tiny fixtures.
- `--soil_intensity`, `--soil_mgf`, `--sample_metadata`, `--enable_rie` and the
  matching/evaluation tuning parameters in `nextflow_schema.json`.
- `docs/usage.md` and `docs/output.md`.

### `Fixed`

- `--input` is no longer required; apply mode starts without a samplesheet.
- Default internal standard is PC 15:0–18:1(d7) (was LPE 18:0(d7)).
- The L3 RIE correction is opt-in (`--enable_rie`, experimental) rather than
  always applied; it is not part of the published method.
- Removed `--skip_gnps`, which no module used.

### `Dependencies`

### `Deprecated`
