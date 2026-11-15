//
// PROCESS — raw files to quality-checked feature tables (stages 1–3).
//
//   sample sheet (file, polarity, role, ...)
//     │
//     ├─ .raw  ──► CONVERT_RAW (stage 1) ──┐
//     └─ .mzML ───────────────────────────┤
//                                          ▼
//                         per polarity: MZMINE_PROCESS (stage 2)
//                                          ▼
//                         per polarity: QC_REPORT (stage 3)
//
// Later stages (annotation, quantification and recovery, report, Atlas) read
// the outputs of this workflow.
//

include { CONVERT_RAW    } from '../modules/local/convert_raw.nf'
include { MZMINE_PROCESS } from '../modules/local/mzmine_process.nf'
include { QC_REPORT      } from '../modules/local/qc_report.nf'


workflow PROCESS {
    take:
    ch_samplesheet // channel: [ meta, file ] from the validated sample sheet

    main:
    ch_versions = channel.empty()

    if (!params.input) {
        error "Process mode needs --input <samplesheet.csv>. See docs/usage.md."
    }

    // --- Stage 1: convert Thermo raw files; mzML files pass through --------
    ch_split = ch_samplesheet.branch { _meta, f ->
        raw:  f.name.toLowerCase().endsWith('.raw')
        mzml: true
    }
    CONVERT_RAW(ch_split.raw)
    ch_versions = ch_versions.mix(CONVERT_RAW.out.versions.first())

    ch_mzml = ch_split.mzml.mix(CONVERT_RAW.out.mzml)

    // --- Stage 2: one MZmine run per polarity -------------------------------
    ch_by_polarity = ch_mzml
        .map { meta, f -> [ meta.polarity, f ] }
        .groupTuple()

    ch_mzmine_in = ch_by_polarity.map { polarity, files ->
        def batch = polarity == 'pos' ? params.mzmine_batch_pos : params.mzmine_batch_neg
        [ polarity, files, file(batch, checkIfExists: true) ]
    }
    MZMINE_PROCESS(ch_mzmine_in)
    ch_versions = ch_versions.mix(MZMINE_PROCESS.out.versions.first())

    // --- Stage 3: quality report per polarity -------------------------------
    QC_REPORT(
        MZMINE_PROCESS.out.features.join(ch_by_polarity),
        file(params.input, checkIfExists: true),
        file(params.is_targets, checkIfExists: true)
    )
    ch_versions = ch_versions.mix(QC_REPORT.out.versions.first())

    emit:
    features    = MZMINE_PROCESS.out.features
    gnps_mgf    = MZMINE_PROCESS.out.gnps_mgf
    sirius_mgf  = MZMINE_PROCESS.out.sirius_mgf
    annotations = MZMINE_PROCESS.out.annotations
    qc_summary  = QC_REPORT.out.summary
    versions    = ch_versions
}
