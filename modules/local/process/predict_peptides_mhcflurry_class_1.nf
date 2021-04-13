// Import generic module functions
include { initOptions; saveFiles } from './functions'

params.options = [:]

//TODO: combine in a subflow --> when needs to be removed
process PREDICT_PEPTIDES_MHCFLURRY_CLASS_1 {
    publishDir "${params.outdir}",
    mode: params.publish_dir_mode,
    saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:'class_1_bindings', publish_id:'class_1_bindings') }
    
    conda (params.enable_conda ? "bioconda::fred2=2.0.6 bioconda::mhcflurry=1.4.3 bioconda::mhcnuggets=2.3.2" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/mulled-v2-689ae0756dd82c61400782baaa8a7a1c2289930d:a9e10ca22d4cbcabf6b54f0fb5d766ea16bb171e-0"
    } else {
        container "quay.io/biocontainers/mulled-v2-689ae0756dd82c61400782baaa8a7a1c2289930d:a9e10ca22d4cbcabf6b54f0fb5d766ea16bb171e-0"
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