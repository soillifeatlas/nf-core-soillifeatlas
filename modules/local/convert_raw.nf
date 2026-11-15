//
// CONVERT_RAW — stage 1: Thermo .raw to indexed, centroided mzML.
//
// ThermoRawFileParser (Hulstaert et al. 2020, J. Proteome Res.) runs on Linux
// and macOS without the vendor library installed. Vendor peak picking stays on
// (the default), so spectra come out centroided, which MZmine expects.
// Set --thermorawfileparser to the executable's path if it is not on PATH.
//
process CONVERT_RAW {
    tag "${meta.id}"
    label 'process_low'

    input:
    tuple val(meta), path(raw)

    output:
    tuple val(meta), path("${meta.id}.mzML"), emit: mzml
    path 'versions.yml',                    emit: versions

    script:
    def trfp = params.thermorawfileparser ?: 'ThermoRawFileParser'
    def args = task.ext.args ?: ''
    """
    ${trfp} -i=${raw} -o=. -f=2 ${args}
    if [ ! -f "${meta.id}.mzML" ]; then
        mv ./*.mzML "${meta.id}.mzML"
    fi
    test -s "${meta.id}.mzML"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        thermorawfileparser: \$(${trfp} --version 2>/dev/null | head -1)
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.mzML
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        thermorawfileparser: stub
    END_VERSIONS
    """
}
