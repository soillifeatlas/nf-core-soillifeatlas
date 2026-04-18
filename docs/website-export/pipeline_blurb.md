# Stage 7 — Automated pipeline

The soillifeatlas pipeline converts the post-mass-spectrometry workflow from
Samrat et al. 2025 (Soil Biol Biochem, doi:10.1016/j.soilbio.2025.109892)
into a reproducible, one-command Nextflow workflow. It reads raw intensity
data and spectral MS2, applies the analysis-19 IS + RIE (floor=0.20) +
Archaea reference-filter correction pipeline, and decomposes soil samples
against the atlas SIMPER fingerprint using four methods in parallel.

- **Reproducibility:** every run asserts a SHA256-pinned atlas and emits a
  `provenance.yaml` capturing git SHA, atlas DOI, container digests, and
  params. Nightly CI proves the pipeline still builds and runs end-to-end.
- **Deployment:** runs on a laptop (Docker), on LISC SLURM (Apptainer), and
  on Google Cloud Batch (stub — post-funding). Pinned nf-core DSL2 template,
  MIT license, GitHub Actions CI.
- **Status:** v0.1-dev ships today. v0.1.0 tag lands when the reference
  atlas gets its Zenodo DOI.

[→ Repository](https://github.com/soillifeatlas/nf-core-soillifeatlas)
[→ Architecture](https://github.com/soillifeatlas/nf-core-soillifeatlas/blob/main/docs/ARCHITECTURE.md)
