process OPENMS_PERCOLATORADAPTER {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::openms-thirdparty=3.1.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/openms-thirdparty:3.1.0--h9ee0642_3' :
        'biocontainers/openms-thirdparty:3.1.0--h9ee0642_3' }"

    input:
        tuple val(meta), path(merged_with_features)

    output:
        tuple val(meta), path("*.idXML")                         , emit: idxml
        tuple val(meta), path("*_percolator_feature_weights.tsv"), emit: feature_weights, optional: true
        tuple val(meta), path("*_percolator_feature_weights.tsv"), emit: feature_weights, optional: true

    when:
        task.ext.when == null || task.ext.when

    script:
    def prefix           = task.ext.prefix ?: "${meta.id}_pout"
    def args             = task.ext.args  ?: ''
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
    def prefix           = task.ext.prefix ?: "${meta.id}_pout"
    def args             = task.ext.args  ?: ''

    """
    touch ${prefix}.idXML

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        PercolatorAdapter: \$(PercolatorAdapter 2>&1 | grep -E '^Version(.*)' | sed 's/Version: //g' | cut -d ' ' -f 1)
        percolator: \$(percolator -h 2>&1 | grep -E '^Percolator version(.*)' | sed 's/Percolator version //g')
    END_VERSIONS
    """
}
