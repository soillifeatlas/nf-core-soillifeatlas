//
// QC_REPORT — stage 3: run-level quality checks for one polarity.
//
// Internal-standard retention time, mass accuracy and peak-area RSD on QC
// injections and standards, pooled-QC feature reproducibility, blank
// contribution, and per-injection summaries in injection order (taken from the
// mzML headers). See lib/soillifeatlas/process_qc.py for the definitions.
//
process QC_REPORT {
    tag "${polarity}"
    label 'process_medium'

    input:
    tuple val(polarity), path(features), path(mzmls, stageAs: 'mzml/*')
    path samplesheet
    path is_targets
    val lib_sig  // libSignature(): package code is part of the task hash

    output:
    tuple val(polarity), path("${prefix}.qc_report.html"),  emit: html
    tuple val(polarity), path("${prefix}.qc_summary.json"), emit: summary
    tuple val(polarity), path("${prefix}.qc_*.tsv"),        emit: tables
    path 'versions.yml',                                     emit: versions

    script:
    prefix = "${params.run_name}_${polarity}"
    def compounds = params.is_compounds ? "--is-compounds '${params.is_compounds}'" : ''
    """
    qc_report.py \\
        --features ${features} \\
        --samplesheet ${samplesheet} \\
        --targets ${is_targets} \\
        --polarity ${polarity} \\
        --mzml mzml/*.mzML \\
        --prefix ${prefix} \\
        --rt-tolerance ${params.qc_rt_tolerance} \\
        --is-rsd ${params.qc_is_rsd} \\
        --mass-accuracy-ppm ${params.qc_mass_accuracy_ppm} \\
        --qc-rsd ${params.qc_feature_rsd} \\
        --qc-pass-fraction ${params.qc_feature_pass_fraction} \\
        --blank-ratio ${params.qc_blank_ratio} \\
        --rt-window ${params.is_rt_window} \\
        ${compounds}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
        soillifeatlas: \$(python3 -c 'import soillifeatlas; print(soillifeatlas.__version__)')
        soillifeatlas_code: ${lib_sig}
    END_VERSIONS
    """

    stub:
    prefix = "${params.run_name}_${polarity}"
    """
    touch ${prefix}.qc_report.html ${prefix}.qc_summary.json ${prefix}.qc_checks.tsv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: stub
    END_VERSIONS
    """
}
