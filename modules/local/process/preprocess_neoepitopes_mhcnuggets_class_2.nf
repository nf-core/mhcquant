// Import generic module functions
include { initOptions; saveFiles } from './functions'

params.options = [:]

def VERSION = '2.3.2'

//TODO: combine in a subflow --> when needs to be removed
process PREPROCESS_NEOEPITOPES_MHCNUGGETS_CLASS_2 {
    conda (params.enable_conda ? "bioconda::mhcnuggets=2.3.2--py_0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/mhcnuggets:2.3.2--py_0"
    } else {
        container "quay.io/biocontainers/mhcnuggets:2.3.2--py_0"
    }

    input:
        tuple val(id), val(Sample), file(neoepitopes)

    output:
        tuple val("$id"), val("$Sample"), file("${Sample}_mhcnuggets_preprocessed"), emit: preprocessed   
        path  "*.version.txt", emit: version

    when:
        params.include_proteins_from_vcf & params.predict_class_2

    script:
        """
            preprocess_neoepitopes_mhcnuggets.py --neoepitopes ${neoepitopes} --output ${Sample}_mhcnuggets_preprocessed

            echo $VERSION > mhcnuggets.version.txt
        """
}
