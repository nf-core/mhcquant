process OPENMS_IDCONFLICTRESOLVER {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/openms:3.4.1--h81ffffe_1' :
        'biocontainers/openms:3.4.1--h81ffffe_1' }"

    input:
    tuple val(meta), path(consensus)

    output:
    tuple val(meta), path("*.consensusXML"), emit: consensusxml
    tuple val("${task.process}"), val('openms'), eval("FileInfo --help 2>&1 | grep -E '^Version' | sed 's/^.*Version: //; s/-.*\$//' | sed 's/ -*//; s/ .*\$//'"), emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}_resolved"

    """
    IDConflictResolver -in $consensus \\
        -out ${prefix}.consensusXML \\
        -threads $task.cpus
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}_resolved"

    """
    touch ${prefix}.consensusXML
    """
}
