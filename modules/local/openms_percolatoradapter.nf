process OPENMS_PERCOLATORADAPTER {
    tag "$meta.id"
    label 'process_high'

    conda (params.enable_conda ? "bioconda::openms-thirdparty=2.6.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/openms-thirdparty:2.6.0--0' :
        'quay.io/biocontainers/openms-thirdparty:2.6.0--0' }"

    input:
        tuple val(meta), path(psm)

    output:
        tuple val(meta), path("*.idXML"), emit: idxml
        path "versions.yml"             , emit: versions

    script:
        def prefix           = task.ext.prefix ?: "${meta.id}"
        def args             = task.ext.args  ?: ''
        def klammer          = (params.description_correct_features > 0 && params.klammer) ? "-klammer" : ""

        """
        OMP_NUM_THREADS=$task.cpus \\
        PercolatorAdapter -in $psm \\
            -out ${prefix}.idXML \\
            $klammer \\
            $args

        cat <<-END_VERSIONS > versions.yml
        ${task.process}:
            openms-thirdparty: \$(echo \$(FileInfo --help 2>&1) | sed 's/^.*Version: //; s/-.*\$//' | sed 's/ -*//; s/ .*\$//')
        END_VERSIONS
        """
}
