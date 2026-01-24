process OPENMS_MAPALIGNERIDENTIFICATION {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/openms:3.5.0--h78fb946_0' :
        'biocontainers/openms:3.5.0--h78fb946_0' }"

    input:
    tuple val(meta), path(idxmls)

    output:
    tuple val(meta), path("*.trafoXML"), emit: trafoxml
    tuple val("${task.process}"), val('openms'), eval("FileInfo --help 2>&1 | sed -nE 's/^Version: ([0-9.]+).*/\\1/p'"), emit: versions_openms, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args      = task.ext.args  ?: ''
    def out_names = idxmls.collect { it.baseName.replace('_fdr_filtered','')+'.trafoXML' }.join(' ')
    
    """
    MapAlignerIdentification -in $idxmls \\
        -trafo_out ${out_names} \\
        $args
    """

    stub:

    """
    touch test1.consensusXML
    touch test2.consensusXML
    """
}
