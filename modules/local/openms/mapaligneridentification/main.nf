process OPENMS_MAPALIGNERIDENTIFICATION {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/openms:3.4.1--h81ffffe_1' :
        'biocontainers/openms:3.4.1--h81ffffe_1' }"

    input:
    tuple val(meta), path(idxmls)

    output:
    tuple val(meta), path("*.trafoXML"), emit: trafoxml
    tuple val("${task.process}"), val('openms'), eval("FileInfo --help 2>&1 | grep -E '^Version' | sed 's/^.*Version: //; s/-.*\$//' | sed 's/ -*//; s/ .*\$//'"), emit: versions, topic: versions

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
