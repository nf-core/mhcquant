// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]

process PREDICT_NEOEPITOPES_MHCFLURRY_CLASS_1 {
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:'class_1_bindings', publish_id:'class_1_bindings') }
    
    conda (params.enable_conda ? "bioconda::mhcflurry=1.4.3--py_0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/mhcflurry:1.4.3--py_0"
    } else {
        container "quay.io/biocontainers/mhcflurry:1.4.3--py_0"
    }


    input:
        tuple val(Sample), val(id), val(allotypes), val(d), path(neoepitopes) 

    output:
        tuple val("$id"), val("$Sample"), path("*_${Sample}_predicted_neoepitopes_class_1.csv"), emit: csv   
        path  "*.version.txt", emit: version

    script:
        """
            mhcflurry-downloads --quiet fetch models_class1
            mhcflurry_neoepitope_binding_prediction.py '${allotypes}' ${neoepitopes} _${Sample}_predicted_neoepitopes_class_1.csv
            mhcflurry-predict --version &> mhcflurry.version.txt
        """
}