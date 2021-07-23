// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

def VERSION = '2.3.2'

process POSTPROCESS_PEPTIDES_MHCNUGGETS_CLASS_2 {
    tag "$meta"
    label 'process_low'

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
        tuple val(meta), path(peptides), path(peptide_to_geneID)

    output:
        tuple val(meta), path('*.csv'), emit: csv   
        path  "*.version.txt", emit: version

    script:
        def prefix = options.suffix ? "${meta.sample}_${options.suffix}" : "${meta.sample}_postprocessed"
    
        """
            postprocess_peptides_mhcnuggets.py --input ${peptides} --peptides_seq_ID ${peptide_to_geneID} --output ${prefix}.csv
            echo $VERSION > mhcnuggets.version.txt
        """
}
