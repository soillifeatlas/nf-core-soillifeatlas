process CORRECT_L5_REF_FILTER {
    tag 'L5_refilter'
    label 'process_low'

    input:
    path simper_atlas
    path archlips_validated

    output:
    path 'simper_atlas_l5_filtered.parquet', emit: filtered
    path 'versions.yml',                     emit: versions

    script:
    def archaea_kingdom = params.archaea_kingdom ?: 'Archaea'
    """
    apply_ref_quality_filter.py \\
        --simper-atlas ${simper_atlas} \\
        --archlips-validated ${archlips_validated} \\
        --archaea-kingdom '${archaea_kingdom}' \\
        --output simper_atlas_l5_filtered.parquet

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
        soillifeatlas: \$(python3 -c 'import soillifeatlas; print(soillifeatlas.__version__)')
    END_VERSIONS
    """
}
