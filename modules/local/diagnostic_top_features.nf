process DIAGNOSTIC_TOP_FEATURES {
    tag 'top_features'
    label 'process_single'

    input:
    path corrected_intensity
    path simper_atlas

    output:
    path 'top_features_per_phylum.tsv', emit: tsv
    path 'versions.yml',                emit: versions

    script:
    def top_n = params.diagnostic_top_n ?: 10
    """
    diagnose_top_features.py \\
        --corrected-intensity ${corrected_intensity} \\
        --simper-atlas ${simper_atlas} \\
        --top-n ${top_n} \\
        --output top_features_per_phylum.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
        soillifeatlas: \$(python3 -c 'import soillifeatlas; print(soillifeatlas.__version__)')
    END_VERSIONS
    """
}
