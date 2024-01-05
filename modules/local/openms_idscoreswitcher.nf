process OPENMS_IDSCORESWITCHER {
    tag "$meta.id"
    label 'process_single'

    conda "bioconda::openms=3.1.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/openms:3.1.0--h8964181_3' :
        'biocontainers/openms:3.1.0--h8964181_3' }"

    input:
        tuple val(meta), path(idxml), path(whitelist)

    output:
        tuple val(meta), path("*.idXML"), path(whitelist), emit: switched_idxml
        path "versions.yml"             , emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        // TODO: fix naming to be more generic
        def prefix           = task.ext.prefix ?: "${meta.id}_${meta.sample}_${meta.condition}_switched"
        def args             = task.ext.args  ?: ''

        """
        IDScoreSwitcher -in $idxml \\
            -out  ${prefix}.idXML \\
            -threads $task.cpus \\
            $args

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            openms: \$(echo \$(FileInfo --help 2>&1) | sed 's/^.*Version: //; s/-.*\$//' | sed 's/ -*//; s/ .*\$//')
        END_VERSIONS
        """
}
