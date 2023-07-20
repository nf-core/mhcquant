process OPENMS_IDFILTER {
    tag "$meta.id"
    label 'process_single'

    conda "bioconda::openms=2.9.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/openms:2.9.1--h135471a_1' :
        'biocontainers/openms:2.9.1--h135471a_1' }"

    input:
        tuple val(meta), path(idxml), file(peptide_filter)

    output:
        tuple val(meta), path("*.idXML"), emit: idxml
        path "versions.yml"             , emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def whitelist        = "$peptide_filter"
        def prefix           = task.ext.prefix ?: "${meta.id}_-_${idxml.baseName}_filtered"
        def args             = task.ext.args  ?: ''

        if (whitelist == "input.2") {
            whitelist = " "
        }

        """
        IDFilter -in $idxml \\
            -out ${prefix}.idXML \\
            -threads $task.cpus \\
            $args \\
            $whitelist

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            openms: \$(echo \$(FileInfo --help 2>&1) | sed 's/^.*Version: //; s/-.*\$//' | sed 's/ -*//; s/ .*\$//')
        END_VERSIONS
        """
}
