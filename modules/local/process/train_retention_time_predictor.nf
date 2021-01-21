// Import generic module functions
include { initOptions; saveFiles } from './functions'

params.options = [:]

process TRAIN_RETENTION_TIME_PREDICTOR {

    input:
    tuple val(id), val(Sample), file(id_files_for_rt_training)

    output:
    tuple val("$id"), val("$Sample"), file("${Sample}_id_files_for_rt_training.txt"), file("${Sample}_id_files_for_rt_training_params.paramXML"), file("${Sample}_id_files_for_rt_training_trainset.txt")

    when:
    params.predict_RT

    script:
    """
    RTModel -in ${id_files_for_rt_training} \\
        -cv:skip_cv \\
        -out ${Sample}_id_files_for_rt_training.txt \\
        -out_oligo_params ${Sample}_id_files_for_rt_training_params.paramXML \\
        -out_oligo_trainset ${Sample}_id_files_for_rt_training_trainset.txt
    """
}