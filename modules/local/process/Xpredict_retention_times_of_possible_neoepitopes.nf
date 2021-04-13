// Import generic module functions
include { initOptions; saveFiles } from './functions'

params.options = [:]

//TODO: combine in a subflow --> when needs to be removed
process PREDICT_RETENTION_TIMES_OF_POSSIBLE_NEOEPITOPES {
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:'RT_prediction', publish_id:'RT_prediction') }

    conda (params.enable_conda ? "bioconda::openms-thirdparty=2.5.0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/openms-thirdparty:2.5.0--6"
    } else {
        container "quay.io/biocontainers/openms-thirdparty:2.5.0--6"
    }

    input:
        tuple val(id), val(Sample), file(txt_file_for_rt_prediction), file(trained_rt_model_file_II), file(trained_rt_param_file_II), file(trained_rt_set_file_II) 

    output:
        tuple val("$Sample"), file("${Sample}_txt_file_for_rt_prediction_RTpredicted.csv"), emit: csv   
        path  "*.version.txt", emit: version

    when:
        params.predict_RT

    script:
    """
        RTPredict -in_text ${txt_file_for_rt_prediction} \\
            -svm_model ${trained_rt_model_file_II} \\
            -in_oligo_params ${trained_rt_param_file_II} \\
            -in_oligo_trainset ${trained_rt_set_file_II} \\
            -out_text:file ${Sample}_txt_file_for_rt_prediction_RTpredicted.csv

        FileInfo --help &> openms.version.txt
    """
}