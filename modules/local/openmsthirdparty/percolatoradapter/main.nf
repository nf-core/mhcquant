process OPENMS_PERCOLATORADAPTER {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/openms-thirdparty:3.4.0--h9ee0642_0' :
        'biocontainers/openms-thirdparty:3.4.0--h9ee0642_0' }"

    input:
    tuple val(meta), path(merged_with_features)

    output:
    tuple val(meta), path("*.idXML")                         , emit: idxml
    tuple val(meta), path("*_percolator_feature_weights.tsv"), emit: feature_weights, optional: true
    path "versions.yml"                                      , emit: versions

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

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        PercolatorAdapter: \$(PercolatorAdapter 2>&1 | grep -E '^Version(.*)' | sed 's/Version: //g' | cut -d ' ' -f 1)
        percolator: \$(percolator -h 2>&1 | grep -E '^Percolator version(.*)' | sed 's/Percolator version //g')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}_pout"

    """
    touch ${prefix}.idXML

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        PercolatorAdapter: \$(PercolatorAdapter 2>&1 | grep -E '^Version(.*)' | sed 's/Version: //g' | cut -d ' ' -f 1)
        percolator: \$(percolator -h 2>&1 | grep -E '^Percolator version(.*)' | sed 's/Percolator version //g')
    END_VERSIONS
    """
}
