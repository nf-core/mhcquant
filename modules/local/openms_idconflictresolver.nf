process OPENMS_IDCONFLICTRESOLVER {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::openms=2.6.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/openms:2.6.0--h4afb90d_0' :
        'quay.io/biocontainers/openms:2.6.0--h4afb90d_0' }"

    input:
        tuple val(meta), path(consensus)

    output:
        tuple val(meta), path("*.consensusXML"), emit: consensusxml
        path "versions.yml"                    , emit: versions

    script:
        def prefix           = task.ext.suffix ? "${meta.id}_${task.ext.suffix}" : "${meta.id}_resolved"

        """
        IDConflictResolver -in $consensus \\
            -out ${prefix}.consensusXML \\
            -threads $task.cpus

        cat <<-END_VERSIONS > versions.yml
        ${task.process}:
            openms: \$(echo \$(FileInfo --help 2>&1) | sed 's/^.*Version: //; s/-.*\$//' | sed 's/ -*//; s/ .*\$//')
        END_VERSIONS
        """
}
