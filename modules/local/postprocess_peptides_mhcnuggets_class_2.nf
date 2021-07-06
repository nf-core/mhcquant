// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]

def VERSION = '2.3.2'

process POSTPROCESS_PEPTIDES_MHCNUGGETS_CLASS_2 {
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:'class_2_bindings', publish_id:'class_2_bindings') }
    
    conda (params.enable_conda ? "bioconda::mhcnuggets=2.3.2--py_0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/mhcnuggets:2.3.2--py_0"
    } else {
        container "quay.io/biocontainers/mhcnuggets:2.3.2--py_0"
    }

    input:
        tuple val(Sample), val(id), path(predicted_peptides), val(d), path(peptide_to_geneID)

    output:
        tuple val("$Sample"), path('*.csv'), emit: csv   
        path  "*.version.txt", emit: version

    script:
        """
            postprocess_peptides_mhcnuggets.py --input ${predicted_peptides} --peptides_seq_ID ${peptide_to_geneID} --output ${Sample}_postprocessed.csv
            echo $VERSION > mhcnuggets.version.txt
        """
}
