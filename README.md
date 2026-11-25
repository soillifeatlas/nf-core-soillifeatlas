# soillifeatlas/nf-core-soillifeatlas

[![CI](https://github.com/soillifeatlas/nf-core-soillifeatlas/actions/workflows/ci.yml/badge.svg)](https://github.com/soillifeatlas/nf-core-soillifeatlas/actions/workflows/ci.yml)
[![Full Demo Reproducibility](https://github.com/soillifeatlas/nf-core-soillifeatlas/actions/workflows/full-demo.yml/badge.svg)](https://github.com/soillifeatlas/nf-core-soillifeatlas/actions/workflows/full-demo.yml)
[![nf-core/tools](https://img.shields.io/badge/nf--core--tools-3.5.2-24B064?logo=nfcore&logoColor=white)](https://nf-co.re/tools)
[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A524.04.0-23aa62.svg)](https://www.nextflow.io/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

<!-- TODO: add Zenodo DOI badge after v0.1.0 release -->

[![Open in GitHub Codespaces](https://img.shields.io/badge/Open_In_GitHub_Codespaces-black?labelColor=grey&logo=github)](https://github.com/codespaces/new/soillifeatlas/nf-core-soillifeatlas)
[![nf-test](https://img.shields.io/badge/unit_tests-nf--test-337ab7.svg)](https://www.nf-test.com)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)
[![Launch on Seqera Platform](https://img.shields.io/badge/Launch%20%F0%9F%9A%80-Seqera%20Platform-%234256e7)](https://cloud.seqera.io/launch?pipeline=https://github.com/soillifeatlas/nf-core-soillifeatlas)

## Current status: v0.1-dev

- **Process mode** (`--mode process`): stages 1–4 (conversion, MZmine, QC report,
  annotation merge) run on real data with `-profile lisc`; reference run: ClimGrass,
  both polarities, commit 2ed5d70 (`paper/`).
- Stage 5 (Q1–Q4 amounts), stage 6 (report) and stage 7 (Atlas matching) exist as
  library code and scripts used for the manuscript (`paper/scripts/`); not yet
  Nextflow processes.

- End-to-end APPLY workflow runs on tiny fixtures (11/11 nf-tests + 32/32 Python tests green)
- CI green on every PR + nightly full-demo on GitHub Actions
- **v0.1.0 release blocked on Zenodo DOI minting for the reference atlas**
- Real-data reproduction of Samrat 2025 kingdom numbers: queued for v0.1.0 nightly CI once atlas DOI is live

See [`docs/ARCHITECTURE.md`](docs/ARCHITECTURE.md) for the 1-page reviewer
artefact (DAG + SOP traceability + deployment matrix + "not vaporware"
contract).

## Quickstart (macOS / Linux, ~20 min)

Check that the pipeline runs on the bundled tiny fixtures (host Python with
`pip install -e .`):

```bash
nextflow run soillifeatlas/nf-core-soillifeatlas -profile test --outdir results
```

Your own atlas-aligned soil data:

```bash
nextflow run soillifeatlas/nf-core-soillifeatlas \
  -profile local \
  --soil_intensity soil_quant.parquet \
  --soil_mgf soil.mgf \
  --sample_metadata sample_metadata.tsv \
  --outdir results
```

> No release tag or Zenodo DOI exists yet, and no module declares a container,
> so `-profile docker` would still use the host Python. See
> [docs/usage.md](docs/usage.md) for inputs and which level of inference your
> data can reach, and [docs/output.md](docs/output.md) for the results.

## Introduction

**soillifeatlas/nf-core-soillifeatlas** applies the Soil Life Atlas, a
cross-kingdom LC–MS/MS lipid reference, to soil lipidomics data. It matches soil
MS/MS spectra to atlas fingerprint features, normalises intensities to the
PC 15:0–18:1(d7) internal standard, and decomposes each sample's lipid signal
into kingdom and phylum shares with four methods, followed by plausibility and
diagnostic checks. The outputs are shares of lipid signal, not biomass.

## Usage

See [docs/usage.md](docs/usage.md).

## Credits

soillifeatlas/soillifeatlas was originally written by Rahul Samrat.

We thank the following people for their extensive assistance in the development of this pipeline:

<!-- TODO nf-core: If applicable, make list of people who have also contributed -->

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

## Citations

<!-- TODO nf-core: Add citation for pipeline after first release. Uncomment lines below and update Zenodo doi and badge at the top of this file. -->
<!-- If you use soillifeatlas/soillifeatlas for your analysis, please cite it using the following doi: [10.5281/zenodo.XXXXXX](https://doi.org/10.5281/zenodo.XXXXXX) -->

<!-- TODO nf-core: Add bibliography of tools and data used in your pipeline -->

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

This pipeline uses code and infrastructure developed and maintained by the [nf-core](https://nf-co.re) community, reused here under the [MIT license](https://github.com/nf-core/tools/blob/main/LICENSE).

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
