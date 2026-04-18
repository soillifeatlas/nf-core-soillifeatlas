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

- End-to-end APPLY workflow runs on tiny fixtures (11/11 nf-tests + 32/32 Python tests green)
- CI green on every PR + nightly full-demo on GitHub Actions
- **v0.1.0 release blocked on Zenodo DOI minting for the reference atlas**
- Real-data reproduction of Samrat 2025 kingdom numbers: queued for v0.1.0 nightly CI once atlas DOI is live

See [`docs/ARCHITECTURE.md`](docs/ARCHITECTURE.md) for the 1-page reviewer
artefact (DAG + SOP traceability + deployment matrix + "not vaporware"
contract).

## Quickstart (macOS / Linux, ~20 min)

Local Docker (once the atlas Zenodo release is minted):

```bash
nextflow run soillifeatlas/nf-core-soillifeatlas \
  -r v0.1.0 \
  -profile docker \
  --mode apply \
  --soil_intensity path/to/soil_quant.parquet \
  --soil_mgf path/to/soil.mgf \
  --sample_metadata path/to/sample_metadata.tsv
```

v0.1 demo (runs today, on tiny fixtures):

```bash
cd nf-core-soillifeatlas
pip install -e .
nf-test test tests/workflows/apply.nf.test --profile local
```

> The `v0.1.0` tag and corresponding Zenodo DOI are not minted yet; the first
> block above is the target invocation for the v0.1.0 release. Use the tiny-
> fixture demo (second block) to exercise the pipeline end-to-end today.

## Introduction

**soillifeatlas/nf-core-soillifeatlas** is a standardized Nextflow pipeline for cross-kingdom soil lipidomics decomposition.

<!-- TODO nf-core:
   Complete this sentence with a 2-3 sentence summary of what types of data the pipeline ingests, a brief overview of the
   major pipeline sections and the types of output it produces. You're giving an overview to someone new
   to nf-core here, in 15-20 seconds. For an example, see https://github.com/nf-core/rnaseq/blob/master/README.md#introduction
-->

<!-- TODO nf-core: Include a figure that guides the user through the major workflow steps. Many nf-core
     workflows use the "tube map" design for that. See https://nf-co.re/docs/guidelines/graphic_design/workflow_diagrams#examples for examples.   -->
<!-- TODO nf-core: Fill in short bullet-pointed list of the default steps in the pipeline -->

## Usage

> [!NOTE]
> If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/usage/installation) on how to set-up Nextflow.

<!-- TODO nf-core: Describe the minimum required steps to execute the pipeline, e.g. how to prepare samplesheets.
     Explain what rows and columns represent. For instance (please edit as appropriate):

First, prepare a samplesheet with your input data that looks as follows:

`samplesheet.csv`:

```csv
sample,fastq_1,fastq_2
CONTROL_REP1,AEG588A1_S1_L002_R1_001.fastq.gz,AEG588A1_S1_L002_R2_001.fastq.gz
```

Each row represents a fastq file (single-end) or a pair of fastq files (paired end).

-->

Now, you can run the pipeline using:

<!-- TODO nf-core: update the following command to include all required parameters for a minimal example -->

```bash
nextflow run soillifeatlas/soillifeatlas \
   -profile <docker/singularity/.../institute> \
   --input samplesheet.csv \
   --outdir <OUTDIR>
```

> [!WARNING]
> Please provide pipeline parameters via the CLI or Nextflow `-params-file` option. Custom config files including those provided by the `-c` Nextflow option can be used to provide any configuration _**except for parameters**_; see [docs](https://nf-co.re/docs/usage/getting_started/configuration#custom-configuration-files).

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
