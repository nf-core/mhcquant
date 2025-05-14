process OPENMS_MAPALIGNERIDENTIFICATION {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/openms:3.4.0--hc77a4c7_0' :
        'biocontainers/openms:3.4.0--hc77a4c7_0' }"

    input:
    tuple val(meta), path(idxmls)

    output:
    tuple val(meta), path("*.trafoXML"), emit: trafoxml
    path "versions.yml"                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args      = task.ext.args  ?: ''
    def out_names = idxmls.collect { it.baseName.replace('_fdr_filtered','')+'.trafoXML' }.join(' ')

    """
    MapAlignerIdentification \\
        -in $idxmls \\
        -trafo_out ${out_names} \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        openms: \$(echo \$(FileInfo --help 2>&1) | sed 's/^.*Version: //; s/-.*\$//' | sed 's/ -*//; s/ .*\$//')
    END_VERSIONS
    """

    stub:

    """
    touch test1.consensusXML
    touch test2.consensusXML

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        openms: \$(echo \$(FileInfo --help 2>&1) | sed 's/^.*Version: //; s/-.*\$//' | sed 's/ -*//; s/ .*\$//')
    END_VERSIONS
    """
}
