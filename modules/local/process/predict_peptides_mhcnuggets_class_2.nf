// Import generic module functions
include { initOptions; saveFiles } from './functions'

params.options = [:]

def VERSION = '2.3.2'

//TODO: combine in a subflow --> when needs to be removed
process PREDICT_PEPTIDES_MHCNUGGETS_CLASS_2 {
    conda (params.enable_conda ? "bioconda::mhcnuggets=2.3.2--py_0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/mhcnuggets:2.3.2--py_0"
    } else {
        container "quay.io/biocontainers/mhcnuggets:2.3.2--py_0"
    }

    input:
        tuple val(Sample), val(id), file(preprocessed_peptides), val(d), val(class_2_alleles)

    output:
        tuple val("$id"), val("$Sample"), file("*_predicted_peptides_class_2"), emit: csv   
        path  "*.version.txt", emit: version

    when:
        params.predict_class_2

    script:
        """
            mhcnuggets_predict_peptides.py --peptides ${preprocessed_peptides} --alleles '${class_2_alleles}' --output _predicted_peptides_class_2

            echo $VERSION > mhcnuggets.version.txt
        """
}