process OPENMS_FEATURELINKERUNLABELEDKD {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/openms-thirdparty:3.4.1--h9ee0642_1' :
        'biocontainers/openms-thirdparty:3.4.1--h9ee0642_1' }"

    input:
    tuple val(meta), path(features)

    output:
    tuple val(meta), path("*.consensusXML"), emit: consensusxml
    tuple val("${task.process}"), val('openms'), eval("FileInfo --help 2>&1 | sed -nE 's/^Version: ([0-9.]+).*/\1/p'"), emit: versions_openms, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}_all_features_merged"

    """
    FeatureLinkerUnlabeledKD -in $features \\
        -out ${prefix}.consensusXML \\
        -threads $task.cpus
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}_all_features_merged"

    """
    touch ${prefix}.consensusXML
    """
}
