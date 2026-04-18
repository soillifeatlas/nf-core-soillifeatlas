process DECOMPOSE {
    tag "${method}"
    label 'process_medium'

    input:
    path soil_intensity
    path atlas_intensity
    path simper_atlas
    path sample_phylum_map
    each method

    output:
    tuple val(method), path("composition_kingdom_${method}.parquet"), emit: kingdom
    tuple val(method), path("composition_phylum_${method}.parquet"),  emit: phylum
    path 'versions.yml',                                              emit: versions

    script:
    def min_phylum_samples = params.min_phylum_samples ?: 2
    """
    run_decomposition.py \\
        --soil-intensity ${soil_intensity} \\
        --atlas-intensity ${atlas_intensity} \\
        --simper-atlas ${simper_atlas} \\
        --sample-phylum-map ${sample_phylum_map} \\
        --method ${method} \\
        --min-phylum-samples ${min_phylum_samples} \\
        --output composition_kingdom_${method}.parquet \\
        --output-phylum composition_phylum_${method}.parquet

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
        soillifeatlas: \$(python3 -c 'import soillifeatlas; print(soillifeatlas.__version__)')
    END_VERSIONS
    """
}
