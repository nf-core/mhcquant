// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process OPENMS_RTPREDICT {
    tag "$meta.id"
    label 'process_low'

    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:'RT_prediction', publish_id:'RT_prediction') }

    conda (params.enable_conda ? "bioconda::openms-thirdparty=2.6.0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/openms-thirdparty:2.6.0--0"
    } else {
        container "quay.io/biocontainers/openms-thirdparty:2.6.0--0"
    }

    input:
        tuple val(meta), path(idxml), path(rt_model), path(rt_params), path(trainset)

    output:
        tuple val(meta), path("*.csv"), emit: csv
        path "versions.yml"           , emit: versions

    script:
        def software = getSoftwareName(task.process)
        def prefix = options.suffix ? "${meta.sample}_${options.suffix}" : "${meta.sample}_RTpredicted"

        """
        RTPredict -in_id $idxml \\
            -svm_model $rt_model \\
            -in_oligo_params $rt_params \\
            -in_oligo_trainset $trainset \\
            -out_text:file ${prefix}.csv
            
        cat <<-END_VERSIONS > versions.yml
        ${getProcessName(task.process)}:
            openms-thirdparty: \$(echo \$(FileInfo --help 2>&1) | sed 's/^.*Version: //; s/-.*\$//' | sed 's/ -*//; s/ .*\$//')
        END_VERSIONS
        """
}
