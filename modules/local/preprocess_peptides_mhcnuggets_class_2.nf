// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]

def VERSION = '2.3.2'

process PREPROCESS_PEPTIDES_MHCNUGGETS_CLASS_2 {
    conda (params.enable_conda ? "bioconda::mhcnuggets=2.3.2--py_0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/mhcnuggets:2.3.2--py_0"
    } else {
        container "quay.io/biocontainers/mhcnuggets:2.3.2--py_0"
    }

    input:
        tuple val(id), val(Sample), path(mztab_file) 

    output:
        tuple val("$id"), val("$Sample"), path("${Sample}_preprocessed_mhcnuggets_peptides"), emit: preprocessed
        tuple val("$id"), val("$Sample"), path('peptide_to_geneID'), emit: geneID   
        path  "*.version.txt", emit: version

    script:
        """
            preprocess_peptides_mhcnuggets.py --mztab ${mztab_file} --output ${Sample}_preprocessed_mhcnuggets_peptides
            echo $VERSION > mhcnuggets.version.txt
        """
}
