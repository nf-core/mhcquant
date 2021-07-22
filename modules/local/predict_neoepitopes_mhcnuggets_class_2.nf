// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

def VERSION = '2.3.2'

process PREDICT_NEOEPITOPES_MHCNUGGETS_CLASS_2 {
    tag "$meta"

    conda (params.enable_conda ? "bioconda::mhcnuggets=2.3.2--py_0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/mhcnuggets:2.3.2--py_0"
    } else {
        container "quay.io/biocontainers/mhcnuggets:2.3.2--py_0"
    }

    input:
        tuple val(meta), path(neoepitopes), val(alleles)

    output:
        tuple val(meta), path("*_predicted_neoepitopes_class_2"), emit: csv   
        path  "*.version.txt", emit: version

    script:
        def prefix = options.suffix ? "${meta}_${options.suffix}" : "${meta}_predicted_neoepitopes_class_2"

        """
            mhcnuggets_predict_peptides.py --peptides ${neoepitopes} --alleles '${alleles}' --output ${prefix}
            echo $VERSION > mhcnuggets.version.txt
        """
}