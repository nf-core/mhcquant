process OPENMS_IDFILTER {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::openms=2.6.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/openms:2.6.0--h4afb90d_0' :
        'quay.io/biocontainers/openms:2.6.0--h4afb90d_0' }"

    input:
        tuple val(meta), path(idxml), file(peptide_filter)

    output:
        tuple val(meta), path("*.idXML"), emit: idxml
        path "versions.yml"             , emit: versions

    script:
        def whitelist        = "$peptide_filter"



        def prefix           = task.ext.suffix ? "${meta.id}_-_${idxml.baseName}_${task.ext.suffix}" : "${meta.id}_-_${idxml.baseName}_filtered"
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
        ${task.process}:
            openms: \$(echo \$(FileInfo --help 2>&1) | sed 's/^.*Version: //; s/-.*\$//' | sed 's/ -*//; s/ .*\$//')
        END_VERSIONS
        """
}
