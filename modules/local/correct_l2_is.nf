process CORRECT_L2_IS {
    tag 'L2_IS'
    label 'process_low'

    input:
    path intensity
    path is_features

    output:
    path 'intensity_l2_is_corrected.parquet', emit: corrected
    path 'versions.yml',                      emit: versions

    script:
    def is_reference   = params.IS_reference_compound ?: 'LPE_18d7'
    def is_spiked_pmol = params.IS_spiked_pmol        ?: 100.0
    """
    apply_is_scaling.py \\
        --intensity ${intensity} \\
        --is-features ${is_features} \\
        --is-reference ${is_reference} \\
        --is-spiked-pmol ${is_spiked_pmol} \\
        --output intensity_l2_is_corrected.parquet

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
        soillifeatlas: \$(python3 -c 'import soillifeatlas; print(soillifeatlas.__version__)')
    END_VERSIONS
    """
}
