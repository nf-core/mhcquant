// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]

//TODO: combine in a subflow --> when needs to be removed
process TRAIN_RETENTION_TIME_PREDICTOR {
    conda (params.enable_conda ? "bioconda::openms-thirdparty=2.5.0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/openms-thirdparty:2.5.0--6"
    } else {
        container "quay.io/biocontainers/openms-thirdparty:2.5.0--6"
    }

    input:
        tuple val(id), val(Sample), file(id_files_for_rt_training)

    output:
        tuple val("$id"), val("$Sample"), file("${Sample}_id_files_for_rt_training.txt"), file("${Sample}_id_files_for_rt_training_params.paramXML"), file("${Sample}_id_files_for_rt_training_trainset.txt"), emit: complete   
        path  "*.version.txt", emit: version

    when:
        params.predict_RT

    script:
    """
        RTModel -in ${id_files_for_rt_training} \\
            -cv:skip_cv \\
            -out ${Sample}_id_files_for_rt_training.txt \\
            -out_oligo_params ${Sample}_id_files_for_rt_training_params.paramXML \\
            -out_oligo_trainset ${Sample}_id_files_for_rt_training_trainset.txt

        FileInfo --help &> openms.version.txt
    """
}