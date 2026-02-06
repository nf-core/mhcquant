process OPENMS_PERCOLATORADAPTER {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/openms-thirdparty:3.5.0--h9ee0642_0' :
        'biocontainers/openms-thirdparty:3.5.0--h9ee0642_0' }"

    input:
    tuple val(meta), path(merged_with_features)

    output:
    tuple val(meta), path("*.idXML")                         , emit: idxml
    tuple val(meta), path("*_percolator_feature_weights.tsv"), emit: feature_weights, optional: true
    tuple val("${task.process}"), val('PercolatorAdapter'), eval("PercolatorAdapter 2>&1 | grep -E '^Version(.*)' | sed 's/Version: //g' | cut -d ' ' -f 1"), topic: versions
    tuple val("${task.process}"), val('percolator'), eval("percolator -h 2>&1 | grep -E '^Percolator version(.*)' | sed 's/Percolator version //g'"), topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args  ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}_pout"

    """
    PercolatorAdapter \\
        -in $merged_with_features \\
        -out ${prefix}.idXML \\
        -threads $task.cpus \\
        $args
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}_pout"

    """
    touch ${prefix}.idXML
    """
}
