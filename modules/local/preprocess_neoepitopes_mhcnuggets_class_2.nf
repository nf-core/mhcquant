// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]

def VERSION = '2.3.2'

process PREPROCESS_NEOEPITOPES_MHCNUGGETS_CLASS_2 {
    conda (params.enable_conda ? "bioconda::mhcnuggets=2.3.2--py_0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/mhcnuggets:2.3.2--py_0"
    } else {
        container "quay.io/biocontainers/mhcnuggets:2.3.2--py_0"
    }

    input:
        tuple val(id), val(Sample), path(neoepitopes)

    output:
        tuple val("$id"), val("$Sample"), path("${Sample}_mhcnuggets_preprocessed"), emit: preprocessed   
        path  "*.version.txt", emit: version

    script:
        """
            preprocess_neoepitopes_mhcnuggets.py --neoepitopes ${neoepitopes} --output ${Sample}_mhcnuggets_preprocessed
            echo $VERSION > mhcnuggets.version.txt
        """
}
