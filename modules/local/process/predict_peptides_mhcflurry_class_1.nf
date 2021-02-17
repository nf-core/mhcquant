// Import generic module functions
include { initOptions; saveFiles } from './functions'

params.options = [:]

//TODO: combine in a subflow --> when needs to be removed
process PREDICT_PEPTIDES_MHCFLURRY_CLASS_1 {
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
        tuple val(Sample), val(id), file(mztab_file), val(d), val(class_1_alleles)

    output:
        tuple val("$Sample"), file("*predicted_peptides_class_1.csv"), emit: csv   
        path  "*.version.txt", emit: version

    when:
        params.predict_class_1

    script:
    """
        mhcflurry-downloads --quiet fetch models_class1
        mhcflurry_predict_mztab.py '${class_1_alleles}' ${mztab_file} ${Sample}_predicted_peptides_class_1.csv

        mhcflurry-predict --version &> mhcflurry.version.txt
    """
}