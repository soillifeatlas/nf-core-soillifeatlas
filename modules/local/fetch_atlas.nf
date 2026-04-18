process FETCH_ATLAS {
    tag "atlas_${params.atlas_version}"
    label 'process_low'

    output:
    path "atlas_fetched", emit: atlas_dir
    path 'versions.yml',  emit: versions

    script:
    def atlas_path_arg = params.atlas_path   ? "--atlas-path ${params.atlas_path}"      : ""
    def atlas_sha_arg  = params.atlas_sha256 ? "--atlas-sha256 ${params.atlas_sha256}"  : ""
    """
    fetch_atlas.py \\
        --atlas-version ${params.atlas_version} \\
        --atlas-doi ${params.atlas_doi} \\
        ${atlas_sha_arg} \\
        --cache-dir \${NXF_CACHE_DIR:-\$HOME/.nextflow/cache}/soillifeatlas_atlas \\
        ${atlas_path_arg} \\
        --output-dir atlas_fetched

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
        soillifeatlas: \$(python3 -c 'import soillifeatlas; print(soillifeatlas.__version__)')
    END_VERSIONS
    """
}
