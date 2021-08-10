// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

def VERSION = '2.3.2'

process PREPROCESS_NEOEPITOPES_MHCNUGGETS_CLASS_2 {
    tag "$meta"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::mhcnuggets=2.3.2--py_0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/mhcnuggets:2.3.2--py_0"
    } else {
        container "quay.io/biocontainers/mhcnuggets:2.3.2--py_0"
    }

    input:
        tuple val(meta), path(neoepitopes)

    output:
        tuple val(meta), path("*${prefix}*"), emit: preprocessed
        path  "*.version.txt", emit: version

    script:
        def prefix = options.suffix ? "${meta}_${options.suffix}" : "${meta}_mhcnuggets_preprocessed"

        """
            preprocess_neoepitopes_mhcnuggets.py --neoepitopes ${neoepitopes} --output ${prefix}
            echo $VERSION > mhcnuggets.version.txt
        """
}
