process DATAMASH_HISTOGRAM {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/datamash:1.1.0--0' :
        'biocontainers/datamash:1.1.0--0' }"

    input:
    tuple val(meta), path(tsv)

    output:
    tuple val(meta), path("*.tsv"), emit: binned_tsv
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    datamash \\
        bin:0.0002 1 \\
        $args \\
        < $tsv \\
    | \\
    datamash \\
        --sort \\
        --group 1 \\
        count 1 \\
        > ${prefix}_binned.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        datamash: \$(datamash --version | grep 'datamash' | cut -d ' ' -f4)
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch ${prefix}_binned.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        datamash: \$(datamash --version | grep 'datamash' | cut -d ' ' -f4)
    END_VERSIONS
    """
}
