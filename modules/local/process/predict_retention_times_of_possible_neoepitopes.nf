// Import generic module functions
include { initOptions; saveFiles } from './functions'

params.options = [:]

process PREDICT_RETENTION_TIMES_OF_POSSIBLE_NEOEPITOPES {
    publishDir "${params.outdir}/RT_prediction/"

    input:
    tuple val(id), val(Sample), file(txt_file_for_rt_prediction), file(trained_rt_model_file_II), file(trained_rt_param_file_II), file(trained_rt_set_file_II) 

    output:
    tuple val("$Sample"), file("${Sample}_txt_file_for_rt_prediction_RTpredicted.csv")

    when:
    params.predict_RT

    script:
    """
    RTPredict -in_text ${txt_file_for_rt_prediction} \\
        -svm_model ${trained_rt_model_file_II} \\
        -in_oligo_params ${trained_rt_param_file_II} \\
        -in_oligo_trainset ${trained_rt_set_file_II} \\
        -out_text:file ${Sample}_txt_file_for_rt_prediction_RTpredicted.csv
    """
}