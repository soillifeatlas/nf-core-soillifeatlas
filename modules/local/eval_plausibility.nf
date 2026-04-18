process EVAL_PLAUSIBILITY {
    tag "${method}"
    label 'process_single'

    input:
    tuple val(method), path(composition_kingdom)
    path expected_ref

    output:
    tuple val(method), path("plausibility_${method}.tsv"), emit: tsv
    path 'versions.yml',                                   emit: versions

    script:
    """
    eval_plausibility.py \\
        --composition-kingdom ${composition_kingdom} \\
        --expected-ref ${expected_ref} \\
        --method ${method} \\
        --output plausibility_${method}.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
        soillifeatlas: \$(python3 -c 'import soillifeatlas; print(soillifeatlas.__version__)')
    END_VERSIONS
    """
}
