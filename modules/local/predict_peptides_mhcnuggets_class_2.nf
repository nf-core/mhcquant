// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]

def VERSIONFRED2 = '2.0.6'
def VERSIONMHCNUGGETS = '2.3.2'

process PREDICT_PEPTIDES_MHCNUGGETS_CLASS_2 {

    conda (params.enable_conda ? "bioconda::fred2=2.0.6 bioconda::mhcflurry=1.4.3 bioconda::mhcnuggets=2.3.2" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/mulled-v2-689ae0756dd82c61400782baaa8a7a1c2289930d:a9e10ca22d4cbcabf6b54f0fb5d766ea16bb171e-0"
    } else {
        container "quay.io/biocontainers/mulled-v2-689ae0756dd82c61400782baaa8a7a1c2289930d:a9e10ca22d4cbcabf6b54f0fb5d766ea16bb171e-0"
    }

    input:
        tuple val(Sample), val(id), path(preprocessed_peptides), val(d), val(class_2_alleles)

    output:
        tuple val("$id"), val("$Sample"), path("*_predicted_peptides_class_2"), emit: csv   
        path  "*.version.txt", emit: version

    script:
        """
            mhcnuggets_predict_peptides.py --peptides ${preprocessed_peptides} --alleles '${class_2_alleles}' --output _predicted_peptides_class_2
            echo $VERSIONFRED2 > fred2.version.txt
            echo $VERSIONMHCNUGGETS > mhcnuggets.version.txt
        """
}