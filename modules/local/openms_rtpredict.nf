process OPENMS_RTPREDICT {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::openms-thirdparty=2.8.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/openms-thirdparty:2.8.0--h9ee0642_2' :
        'quay.io/biocontainers/openms-thirdparty:2.8.0--h9ee0642_2' }"

    input:
        tuple val(meta), path(idxml), path(rt_model), path(rt_params), path(trainset)

    output:
        tuple val(meta), path("*.csv"), emit: csv
        path "versions.yml"           , emit: versions

    script:
        def prefix           = task.ext.prefix ?: "${meta.sample}_RTpredicted"

        """
        RTPredict -in_id $idxml \\
            -svm_model $rt_model \\
            -in_oligo_params $rt_params \\
            -in_oligo_trainset $trainset \\
            -out_text:file ${prefix}.csv

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            openms-thirdparty: \$(echo \$(FileInfo --help 2>&1) | sed 's/^.*Version: //; s/-.*\$//' | sed 's/ -*//; s/ .*\$//')
        END_VERSIONS
        """
}
