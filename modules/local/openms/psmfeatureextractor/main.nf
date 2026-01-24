process OPENMS_PSMFEATUREEXTRACTOR {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/openms:3.5.0--h78fb946_0' :
        'biocontainers/openms:3.5.0--h78fb946_0' }"

    input:
    tuple val(meta), path(idxml), path(feature_file)

    output:
    tuple val(meta), path("*.idXML"), emit: idxml
    tuple val("${task.process}"), val('openms'), eval("FileInfo --help 2>&1 | sed -nE 's/^Version: ([0-9.]+).*/\\1/p'"), emit: versions_openms, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix         = task.ext.prefix ?: "${meta.id}_psm"
    def args           = task.ext.args ?: ''
    def extra_features = ""

    """
    extra_features=\$(awk 'NR > 1 && \$1 !~ /psm_file/ {printf \"%s \", \$2}' ${feature_file})
    
    PSMFeatureExtractor -in $idxml \\
        -out ${prefix}.idXML \\
        -threads $task.cpus \\
        -extra \$extra_features \\
        $args
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}_psm"

    """
    touch ${prefix}.idXML
    """
}
