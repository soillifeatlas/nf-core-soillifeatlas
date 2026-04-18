process SIMPER_MATCH {
    tag "${meta.id}"
    label 'process_medium'

    input:
    tuple val(meta), path(soil_mgf)
    path atlas_mgf
    path simper_atlas

    output:
    tuple val(meta), path("${meta.id}_verified_matches.parquet"), emit: matches
    path 'versions.yml', emit: versions

    script:
    def precursor_ppm     = params.simper_precursor_ppm     ?: 10
    def fragment_tol      = params.simper_fragment_tol      ?: 0.02
    def min_cos           = params.simper_min_cos           ?: 0.7
    def min_matched_peaks = params.simper_min_matched_peaks ?: 4
    """
    simper_match.py \\
        --soil-mgf ${soil_mgf} \\
        --atlas-mgf ${atlas_mgf} \\
        --simper-atlas ${simper_atlas} \\
        --precursor-ppm ${precursor_ppm} \\
        --fragment-tol ${fragment_tol} \\
        --min-cos ${min_cos} \\
        --min-matched-peaks ${min_matched_peaks} \\
        --output-matches ${meta.id}_verified_matches.parquet

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
        soillifeatlas: \$(python3 -c 'import soillifeatlas; print(soillifeatlas.__version__)')
        matchms: \$(python3 -c 'import matchms; print(matchms.__version__)')
    END_VERSIONS
    """
}
