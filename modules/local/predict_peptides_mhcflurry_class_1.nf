// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]

def VERSIONFRED2 = '2.0.6'
def VERSIONMHCNUGGETS = '2.3.2'

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
        tuple val(Sample), val(id), path(mztab_file), val(d), val(class_1_alleles)

    output:
        tuple val("$Sample"), path("*predicted_peptides_class_1.csv"), emit: csv   
        path  "*.version.txt", emit: version

    script:
        """
            mhcflurry-downloads --quiet fetch models_class1
            mhcflurry_predict_mztab.py '${class_1_alleles}' ${mztab_file} ${Sample}_predicted_peptides_class_1.csv
            echo $VERSIONFRED2 > fred2.version.txt
            echo $VERSIONMHCNUGGETS > mhcnuggets.version.txt
            mhcflurry-predict --version &> mhcflurry.version.txt
        """
}