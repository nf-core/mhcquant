process PYOPENMS_CHROMATOGRAMEXTRACTOR {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pyopenms:3.3.0--py313h9b5bd11_0' :
        'biocontainers/pyopenms:3.3.0--py313h9b5bd11_0' }"

    input:
    tuple val(meta), path(mzml)

    output:
    tuple val(meta), path("*.csv")  , emit: csv
    path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args  ?: ''
    def prefix = task.ext.prefix ?: "${mzml.baseName}"

    """
    chromatogram_extractor.py \\
        -in $mzml \\
        -out ${prefix}_chrom.csv \\

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pyOpenMS: \$(pip show pyopenms | grep Version | cut -d ' ' -f 2)
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${mzml.baseName}"

    """
    touch ${prefix}_chrom.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pyOpenMS: \$(pip show pyopenms | grep Version | cut -d ' ' -f 2)
    END_VERSIONS
    """
}
