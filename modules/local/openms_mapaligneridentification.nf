process OPENMS_MAPALIGNERIDENTIFICATION {
    tag "$merge_id.id"
    label 'process_single'

    conda "bioconda::openms=2.9.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/openms:2.9.1--h135471a_1' :
        'biocontainers/openms:2.9.1--h135471a_1' }"

    input:
        tuple val(merge_id), val(meta), path(idxmls)

    output:
        tuple val(meta), path("*.trafoXML"), emit: trafoxml
        path "versions.yml"                , emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def out_names        = idxmls.collect { it.baseName+'.trafoXML' }.join(' ')
        def args             = task.ext.args  ?: ''

        """
        MapAlignerIdentification -in $idxmls \\
            -trafo_out ${out_names} \\
            -debug 1000 \\
            $args

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            openms: \$(echo \$(FileInfo --help 2>&1) | sed 's/^.*Version: //; s/-.*\$//' | sed 's/ -*//; s/ .*\$//')
        END_VERSIONS
        """
}
