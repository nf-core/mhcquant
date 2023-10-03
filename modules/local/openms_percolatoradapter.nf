process OPENMS_PERCOLATORADAPTER {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::openms-thirdparty=3.0.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/openms-thirdparty:3.0.0--h9ee0642_1' :
        'biocontainers/openms-thirdparty:3.0.0--h9ee0642_1' }"

    input:
        tuple val(meta), path(merged_with_features)

    output:
        tuple val(meta), path("*.idXML"), emit: idxml
        path "versions.yml"             , emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def prefix           = task.ext.prefix ?: "${meta.id}_pout"
        def args             = task.ext.args  ?: ''
        def klammer          = (params.description_correct_features > 0 && params.klammer) ? "-klammer" : ""

        """
        OMP_NUM_THREADS=$task.cpus \\
        PercolatorAdapter -in $merged_with_features \\
            -out ${prefix}.idXML \\
            $klammer \\
            $args

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            openms-thirdparty: \$(echo \$(FileInfo --help 2>&1) | sed 's/^.*Version: //; s/-.*\$//' | sed 's/ -*//; s/ .*\$//')
        END_VERSIONS
        """
}
