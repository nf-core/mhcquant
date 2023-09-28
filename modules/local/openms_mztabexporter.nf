process OPENMS_MZTABEXPORTER {
    tag "$meta.id"
    label 'process_single'

    conda "bioconda::openms=3.0.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/openms:3.0.0--h8964181_1' :
        'biocontainers/openms:3.0.0--h8964181_1' }"

    input:
        tuple val(meta), path(in_file)

    output:
        tuple val(meta), path("*.mzTab"), emit: mztab
        path "versions.yml"             , emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def prefix           = task.ext.prefix ?: "${meta.id}"
        def args             = task.ext.args  ?: ''

        """
        MzTabExporter -in $in_file \\
            -out ${prefix}.mzTab \\
            -threads $task.cpus \\
            $args

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            openms: \$(echo \$(FileInfo --help 2>&1) | sed 's/^.*Version: //; s/-.*\$//' | sed 's/ -*//; s/ .*\$//')
        END_VERSIONS
        """
}
