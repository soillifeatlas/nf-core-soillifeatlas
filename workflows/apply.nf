//
// APPLY — v0.1 entry workflow for apply-mode runs.
//
// Fetches the reference atlas bundle (via FETCH_ATLAS, which resolves from
// --atlas-path / cache / Zenodo) and wires its artefacts alongside the user's
// soil MS2 + intensity inputs into the DECOMPOSE_APPLY subworkflow.
//
// The atlas_dir emitted by FETCH_ATLAS is a *directory* channel; this
// workflow is responsible for unpacking the canonical filenames out of it
// and turning them into per-artefact path channels for the subworkflow
// take: block. That unpacking happens here, not in DECOMPOSE_APPLY, so the
// subworkflow stays independently nf-testable on raw fixtures.
//

include { FETCH_ATLAS     } from '../modules/local/fetch_atlas.nf'
include { DECOMPOSE_APPLY } from '../subworkflows/local/decompose_apply.nf'


workflow APPLY {
    main:
    ch_versions = Channel.empty()

    // --- Step 1: Fetch atlas artefacts ---------------------------------------
    FETCH_ATLAS()
    ch_versions = ch_versions.mix(FETCH_ATLAS.out.versions)

    // --- Step 2: Resolve soil-side inputs ------------------------------------
    // v0.1 uses explicit --soil_intensity / --soil_mgf / --sample_metadata
    // flags rather than a samplesheet; a samplesheet driver lands in v0.2.
    def soil_intensity  = file(params.soil_intensity,  checkIfExists: true)
    def soil_mgf        = file(params.soil_mgf,        checkIfExists: true)
    def sample_metadata = file(params.sample_metadata, checkIfExists: true)

    // --- Step 3: Unpack atlas_dir into per-artefact channels -----------------
    // FETCH_ATLAS.out.atlas_dir is a queue channel with a single directory
    // element. Use `.first()` to promote it to a value channel so each
    // downstream `.map { d -> file("${d}/<name>") }` can consume it.
    //
    // Artefact names below match the canonical atlas bundle layout (see
    // scripts/build_tiny_atlas_bundle.py for the full list).
    ch_atlas = FETCH_ATLAS.out.atlas_dir.first()

    ch_atlas_mgf            = ch_atlas.map { d -> file("${d}/atlas_consensus.mgf",               checkIfExists: true) }
    ch_atlas_intensity      = ch_atlas.map { d -> file("${d}/consensus_aligned_table.parquet",   checkIfExists: true) }
    ch_simper_atlas         = ch_atlas.map { d -> file("${d}/simper_fingerprint_atlas.parquet",  checkIfExists: true) }
    ch_archlips_validated   = ch_atlas.map { d -> file("${d}/archlips_validated_features.csv",   checkIfExists: true) }
    ch_is_features          = ch_atlas.map { d -> file("${d}/equisplash_IS_masses_POS.csv",      checkIfExists: true) }
    ch_rie_table            = ch_atlas.map { d -> file("${d}/rie_table_s10.csv",                 checkIfExists: true) }
    ch_annotation           = ch_atlas.map { d -> file("${d}/consensus_unified_annotations.csv", checkIfExists: true) }
    ch_sample_phylum_map    = ch_atlas.map { d -> file("${d}/sample_phylum_map.csv",             checkIfExists: true) }
    ch_expected_kingdom_ref = ch_atlas.map { d -> file("${d}/expected_kingdom_composition.csv",  checkIfExists: true) }

    // --- Step 4: Run decomposition subworkflow -------------------------------
    DECOMPOSE_APPLY(
        Channel.value(soil_mgf),
        Channel.value(soil_intensity),
        ch_atlas_mgf,
        ch_atlas_intensity,
        ch_simper_atlas,
        ch_archlips_validated,
        ch_is_features,
        ch_rie_table,
        ch_annotation,
        ch_sample_phylum_map,
        ch_expected_kingdom_ref,
        Channel.value(sample_metadata),
    )
    ch_versions = ch_versions.mix(DECOMPOSE_APPLY.out.versions)

    emit:
    compositions_kingdom = DECOMPOSE_APPLY.out.compositions_kingdom
    compositions_phylum  = DECOMPOSE_APPLY.out.compositions_phylum
    plausibility         = DECOMPOSE_APPLY.out.plausibility
    treatment            = DECOMPOSE_APPLY.out.treatment
    diagnostic           = DECOMPOSE_APPLY.out.diagnostic
    versions             = ch_versions
}
