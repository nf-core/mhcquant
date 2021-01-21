// Import generic module functions
include { initOptions; saveFiles } from './functions'

params.options = [:]

process PREDICT_RETENTION_TIMES_OF_FOUND_PEPTIDES {
    publishDir "${params.outdir}/RT_prediction/"

    input:
    tuple val(id), val(Sample), file(id_files_for_rt_prediction), file(trained_rt_model_file), file(trained_rt_param_file), file(trained_rt_set_file)

    output:
    tuple val("$Sample"), file("${Sample}_id_files_for_rt_prediction_RTpredicted.csv")

    when:
    params.predict_RT

    script:
    """
    RTPredict -in_id ${id_files_for_rt_prediction} \\
        -svm_model ${trained_rt_model_file} \\
        -in_oligo_params ${trained_rt_param_file} \\
        -in_oligo_trainset ${trained_rt_set_file} \\
        -out_text:file ${Sample}_id_files_for_rt_prediction_RTpredicted.csv
    """
}