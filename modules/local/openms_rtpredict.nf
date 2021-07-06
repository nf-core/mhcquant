// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options    = initOptions(params.options)

process OPENMS_RTPREDICT {
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:'RT_prediction', publish_id:'RT_prediction') }

    conda (params.enable_conda ? "bioconda::openms=2.5.0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/openms:2.5.0--h4afb90d_6"
    } else {
        container "quay.io/biocontainers/openms:2.5.0--h4afb90d_6"
    }

    input:
        tuple val(id), val(Sample), path(file_rt_prediction), path(trained_rt_model), path(trained_rt_param), path(trained_rt_set)

    output:
        tuple val("$Sample"), path("${Sample}_id_files_for_rt_prediction_RTpredicted.csv"), emit: csv   
        path  "*.version.txt", emit: version

    script:
        def software = getSoftwareName(task.process)
        def prefix = options.suffix ? "${Sample}_${options.suffix}" : "${Sample}_file_for_rt_prediction_RTpredicted"

        """
            RTPredict -in_id ${file_rt_prediction} \\
                -svm_model ${trained_rt_model} \\
                -in_oligo_params ${trained_rt_param} \\
                -in_oligo_trainset ${trained_rt_set} \\
                -out_text:file ${Sample}_id_files_for_rt_prediction_RTpredicted.csv
            echo \$(FileInfo --help 2>&1) | sed 's/^.*Version: //; s/ .*\$//' &> ${software}.version.txt
        """
}