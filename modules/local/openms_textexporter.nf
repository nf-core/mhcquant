process OPENMS_TEXTEXPORTER {
    tag "$meta.id"
    label 'process_single'

    conda "bioconda::openms=2.9.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/openms:2.9.1--h135471a_1' :
        'quay.io/biocontainers/openms:2.9.1--h135471a_1' }"

    input:
        tuple val(meta), path(consensus_resolved)

    output:
        tuple val(meta), path("*.tsv"), emit: tsv
        path "versions.yml"           , emit: versions

    when:
        task.ext.when == null || task.ext.when

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
