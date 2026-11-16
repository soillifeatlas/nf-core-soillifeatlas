process EVAL_TREATMENT {
    tag 'treatment'
    label 'process_single'

    input:
    path composition_phylum
    path sample_metadata

    output:
    path 'treatment_effects.tsv', emit: tsv
    path 'versions.yml',          emit: versions

    script:
    def contrast = params.treatment_contrast ?: 'drought'
    def contrast_list = (contrast instanceof List ? contrast : contrast.tokenize(',')*.trim()).join(' ')
    """
    eval_treatment_effects.py \\
        --composition-phylum ${composition_phylum} \\
        --sample-metadata ${sample_metadata} \\
        --contrast ${contrast_list} \\
        --output treatment_effects.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
        soillifeatlas: \$(python3 -c 'import soillifeatlas; print(soillifeatlas.__version__)')
    END_VERSIONS
    """
}
