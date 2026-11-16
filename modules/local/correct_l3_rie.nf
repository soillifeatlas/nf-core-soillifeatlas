process CORRECT_L3_RIE {
    tag 'L3_RIE'
    label 'process_low'

    input:
    path intensity
    path rie_table
    path annotation

    output:
    path 'intensity_l3_rie_corrected.parquet', emit: corrected
    path 'versions.yml',                       emit: versions

    script:
    def rie_floor    = params.RIE_floor    ?: 0.20
    def rie_ceiling  = params.RIE_ceiling  ?: 100.0
    def fallback_rie = params.fallback_RIE ?: 1.0
    """
    apply_rie_correction.py \\
        --intensity ${intensity} \\
        --rie-table ${rie_table} \\
        --annotation ${annotation} \\
        --rie-floor ${rie_floor} \\
        --rie-ceiling ${rie_ceiling} \\
        --fallback-rie ${fallback_rie} \\
        --output intensity_l3_rie_corrected.parquet

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
        soillifeatlas: \$(python3 -c 'import soillifeatlas; print(soillifeatlas.__version__)')
    END_VERSIONS
    """
}
