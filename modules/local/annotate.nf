//
// ANNOTATE — stage 4: merge lipid annotations into one table per polarity.
//
// Always uses the MZmine lipid annotations in the feature table. Plug-ins are
// optional files supplied by the user: a LipidSearch export, an MS-DIAL
// alignment export, spectral library matches (GNPS results format) and a
// SIRIUS/CANOPUS summary. One written priority rule assigns Gold, Silver or
// Bronze and the quantification route (Q1-Q3). See
// lib/soillifeatlas/annotate.py for the definitions.
//
process ANNOTATE {
    tag "${polarity}"
    label 'process_low'

    input:
    tuple val(polarity), path(features)
    path is_targets
    path lipidsearch, stageAs: 'plugins/lipidsearch/*'
    path msdial,      stageAs: 'plugins/msdial/*'
    path library,     stageAs: 'plugins/library/*'
    path sirius,      stageAs: 'plugins/sirius/*'
    val lib_sig  // libSignature(): package code is part of the task hash

    output:
    tuple val(polarity), path("${prefix}.annotations.tsv"),        emit: annotations
    tuple val(polarity), path("${prefix}.annotation_sources.tsv"), emit: sources
    tuple val(polarity), path("${prefix}.annotation_summary.json"), emit: summary
    path 'versions.yml',                                           emit: versions

    script:
    prefix = "${params.run_name}_${polarity}"
    // stageAs wraps the path in a list, so test the staged path text, not .name
    def plug = { f, flag -> f.toString().contains('NO_FILE') ? '' : "${flag} ${f}" }
    """
    annotate.py \\
        --features ${features} \\
        --polarity ${polarity} \\
        --prefix ${prefix} \\
        --is-targets ${is_targets} \\
        --mz-ppm ${params.annotate_mz_ppm} \\
        --rt-window ${params.annotate_rt_window} \\
        --library-min-score ${params.annotate_library_min_score} \\
        --library-min-peaks ${params.annotate_library_min_peaks} \\
        ${plug(lipidsearch, '--lipidsearch')} \\
        ${plug(msdial, '--msdial')} \\
        ${plug(library, '--library')} \\
        ${plug(sirius, '--sirius')}

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
    touch ${prefix}.annotations.tsv ${prefix}.annotation_sources.tsv ${prefix}.annotation_summary.json
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: stub
    END_VERSIONS
    """
}
