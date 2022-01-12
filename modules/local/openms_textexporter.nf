process OPENMS_TEXTEXPORTER {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::openms=2.6.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/openms:2.6.0--h4afb90d_0' :
        'quay.io/biocontainers/openms:2.6.0--h4afb90d_0' }"

    input:
        tuple val(meta), path(consensus_resolved)

    output:
        tuple val(meta), path("*.tsv"), emit: tsv
        path "versions.yml"           , emit: versions

    script:
        def prefix           = task.ext.prefix ?: "${meta.id}"
        def args             = task.ext.args  ?: ''

        """
        TextExporter -in $consensus_resolved \\
            -out ${prefix}.tsv \\
            -threads $task.cpus \\
            -id:add_hit_metavalues 0 \\
            -id:peptides_only \\
            $args

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            openms: \$(echo \$(FileInfo --help 2>&1) | sed 's/^.*Version: //; s/-.*\$//' | sed 's/ -*//; s/ .*\$//')
        END_VERSIONS
        """
}
