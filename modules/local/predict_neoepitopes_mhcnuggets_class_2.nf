// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]

def VERSION = '2.3.2'

process PREDICT_NEOEPITOPES_MHCNUGGETS_CLASS_2 {
    conda (params.enable_conda ? "bioconda::mhcnuggets=2.3.2--py_0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/mhcnuggets:2.3.2--py_0"
    } else {
        container "quay.io/biocontainers/mhcnuggets:2.3.2--py_0"
    }

    input:
        tuple val(Sample), val(id), path(preprocessed_neoepitopes), val(d), val(cl_2_alleles)

    output:
        tuple val("$id"), val("$Sample"), path("*_predicted_neoepitopes_class_2"), emit: csv   
        path  "*.version.txt", emit: version

    script:
    """
        mhcnuggets_predict_peptides.py --peptides ${preprocessed_neoepitopes} --alleles '${cl_2_alleles}' --output _predicted_neoepitopes_class_2
        echo $VERSION > mhcnuggets.version.txt
    """
}