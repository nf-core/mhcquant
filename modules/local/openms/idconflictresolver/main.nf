process OPENMS_IDCONFLICTRESOLVER {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/openms:3.4.0--hc77a4c7_0' :
        'biocontainers/openms:3.4.0--hc77a4c7_0' }"

    input:
    tuple val(meta), path(consensus)

    output:
    tuple val(meta), path("*.consensusXML"), emit: consensusxml
    path "versions.yml"                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}_resolved"

    """
    IDConflictResolver \\
        -in $consensus \\
        -out ${prefix}.consensusXML \\
        -threads $task.cpus

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        openms: \$(echo \$(FileInfo --help 2>&1) | sed 's/^.*Version: //; s/-.*\$//' | sed 's/ -*//; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}_resolved"

    """
    touch ${prefix}.consensusXML

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        openms: \$(echo \$(FileInfo --help 2>&1) | sed 's/^.*Version: //; s/-.*\$//' | sed 's/ -*//; s/ .*\$//')
    END_VERSIONS
    """
}
