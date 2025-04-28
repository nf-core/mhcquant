process OPENMS_PSMFEATUREEXTRACTOR {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::openms=3.3.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/openms:3.3.0--h0656172_8' :
        'biocontainers/openms:3.3.0--h0656172_8' }"

    input:
    tuple val(meta), path(idxml), path(feature_file)

    output:
    tuple val(meta), path("*.idXML"), emit: idxml
    path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix         = task.ext.prefix ?: "${meta.id}_psm"
    def args           = task.ext.args ?: ''
    def extra_features = ""

    """
    extra_features=\$(awk 'NR > 1 && \$1 !~ /psm_file/ {printf \"%s \", \$2}' ${feature_file})

    PSMFeatureExtractor \\
        -in $idxml \\
        -out ${prefix}.idXML \\
        -threads $task.cpus \\
        -extra \$extra_features \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        openms: \$(echo \$(FileInfo --help 2>&1) | sed 's/^.*Version: //; s/-.*\$//' | sed 's/ -*//; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}_psm"

    """
    touch ${prefix}.idXML

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        openms: \$(echo \$(FileInfo --help 2>&1) | sed 's/^.*Version: //; s/-.*\$//' | sed 's/ -*//; s/ .*\$//')
    END_VERSIONS
    """
}
