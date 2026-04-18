//
// DECOMPOSE_APPLY subworkflow — wires the 8 Phase-3 modules into the v0.1
// apply-mode pipeline: MS2 match -> corrections -> decomposition (x4 methods)
// -> evaluation + diagnostics.
//
// For v0.1 the subworkflow takes explicit file channels per input. The
// pipeline entry point (main.nf / workflows/soillifeatlas.nf) is responsible
// for unpacking the FETCH_ATLAS output directory into these channels — that
// indirection belongs there, not here, so the subworkflow stays easy to
// nf-test on the tiny fixtures.
//

include { SIMPER_MATCH            } from '../../modules/local/simper_match.nf'
include { CORRECT_L2_IS           } from '../../modules/local/correct_l2_is.nf'
include { CORRECT_L3_RIE          } from '../../modules/local/correct_l3_rie.nf'
include { CORRECT_L5_REF_FILTER   } from '../../modules/local/correct_l5_refilter.nf'
include { DECOMPOSE               } from '../../modules/local/decompose.nf'
include { EVAL_PLAUSIBILITY       } from '../../modules/local/eval_plausibility.nf'
include { EVAL_TREATMENT          } from '../../modules/local/eval_treatment.nf'
include { DIAGNOSTIC_TOP_FEATURES } from '../../modules/local/diagnostic_top_features.nf'


workflow DECOMPOSE_APPLY {
    take:
    soil_mgf                // path  — soil MS2 MGF
    soil_intensity          // path  — soil quant parquet
    atlas_mgf               // path  — atlas consensus MS2 MGF
    atlas_intensity         // path  — atlas consensus intensity parquet
    simper_atlas            // path  — SIMPER fingerprint atlas parquet
    archlips_validated      // path  — ArchLips-validated feature CSV
    is_features             // path  — IS feature CSV (for L2 scaling)
    rie_table               // path  — RIE lookup table CSV
    annotation              // path  — feature -> class/adduct annotation CSV
    sample_phylum_map       // path  — atlas sample_id -> phylum CSV
    expected_kingdom_ref    // path  — expected kingdom composition CSV
    sample_metadata         // path  — soil sample metadata TSV (treatment assignments)

    main:
    ch_versions = Channel.empty()

    // --- Step 1: Layer 5 filter on atlas SIMPER ----------------------------
    CORRECT_L5_REF_FILTER(simper_atlas, archlips_validated)
    ch_versions = ch_versions.mix(CORRECT_L5_REF_FILTER.out.versions)

    // --- Step 2: MS2 cosine match soil vs atlas, joined against SIMPER ------
    ch_soil_mgf = soil_mgf.map { f -> [[id: 'demo'], f] }
    SIMPER_MATCH(ch_soil_mgf, atlas_mgf, CORRECT_L5_REF_FILTER.out.filtered)
    ch_versions = ch_versions.mix(SIMPER_MATCH.out.versions)

    // --- Step 3: L2 IS normalization on soil intensity ----------------------
    CORRECT_L2_IS(soil_intensity, is_features)
    ch_versions = ch_versions.mix(CORRECT_L2_IS.out.versions)

    // --- Step 4: L3 RIE correction -----------------------------------------
    CORRECT_L3_RIE(CORRECT_L2_IS.out.corrected, rie_table, annotation)
    ch_versions = ch_versions.mix(CORRECT_L3_RIE.out.versions)

    // --- Step 5: Decompose via all 4 methods in parallel (`each`) ----------
    DECOMPOSE(
        CORRECT_L3_RIE.out.corrected,
        atlas_intensity,
        CORRECT_L5_REF_FILTER.out.filtered,
        sample_phylum_map,
        params.decomposition_methods,
    )
    ch_versions = ch_versions.mix(DECOMPOSE.out.versions)

    // --- Step 6: Plausibility per method -----------------------------------
    // `expected_kingdom_ref` is a single-file queue channel, but we need the
    // same reference delivered alongside each of the 4 per-method kingdom
    // compositions. Convert to a value channel via `.first()` so Nextflow
    // broadcasts it for every iteration of DECOMPOSE.out.kingdom.
    EVAL_PLAUSIBILITY(DECOMPOSE.out.kingdom, expected_kingdom_ref.first())
    ch_versions = ch_versions.mix(EVAL_PLAUSIBILITY.out.versions)

    // --- Step 7: Treatment effects on the primary method's phylum comp -----
    ch_primary_phylum = DECOMPOSE.out.phylum
        .filter { method, _file -> method == params.primary_method }
        .map { _method, f -> f }
    EVAL_TREATMENT(ch_primary_phylum, sample_metadata)
    ch_versions = ch_versions.mix(EVAL_TREATMENT.out.versions)

    // --- Step 8: Top-features diagnostic -----------------------------------
    DIAGNOSTIC_TOP_FEATURES(CORRECT_L3_RIE.out.corrected, CORRECT_L5_REF_FILTER.out.filtered)
    ch_versions = ch_versions.mix(DIAGNOSTIC_TOP_FEATURES.out.versions)

    emit:
    matches              = SIMPER_MATCH.out.matches
    intensity_l2         = CORRECT_L2_IS.out.corrected
    intensity_l3         = CORRECT_L3_RIE.out.corrected
    simper_filtered      = CORRECT_L5_REF_FILTER.out.filtered
    compositions_kingdom = DECOMPOSE.out.kingdom
    compositions_phylum  = DECOMPOSE.out.phylum
    plausibility         = EVAL_PLAUSIBILITY.out.tsv
    treatment            = EVAL_TREATMENT.out.tsv
    diagnostic           = DIAGNOSTIC_TOP_FEATURES.out.tsv
    versions             = ch_versions
}
