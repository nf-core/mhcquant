process OPENMS_IDRIPPER {
    tag "${meta.id}"
    label 'process_single'

    conda "bioconda::openms=3.0.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/openms:3.0.0--h8964181_1' :
        'biocontainers/openms:3.0.0--h8964181_1' }"

    input:
        tuple val(meta), path(merged_idxml)

    output:
        tuple val(meta), path("*.idXML"), emit: ripped
        path "versions.yml"             , emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def args             = task.ext.args  ?: ''

        """
        IDRipper -in $merged_idxml \\
            -out . \\
            -threads $task.cpus \\
            $args

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            openms: \$(echo \$(FileInfo --help 2>&1) | sed 's/^.*Version: //; s/-.*\$//' | sed 's/ -*//; s/ .*\$//')
        END_VERSIONS
        """
}
