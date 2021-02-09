// Import generic module functions
include { initOptions; saveFiles } from './functions'

params.options = [:]

process INDEX_PEPTIDES {
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:'') }

    conda (params.enable_conda ? "bioconda::openms-thirdparty=2.5.0" : null)
    
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/openms-thirdparty:2.5.0--6"
    } else {
        container "quay.io/biocontainers/openms-thirdparty:2.5.0--6"
    }

    input:
        tuple val(Sample), val(id), val(Condition), file(id_file), val(d), file(fasta_decoy)

    output:
        tuple val("$id"), val("$Sample"), val("$Condition"), file("${Sample}_${Condition}_${id}_idx.idXML")

    script:
    """
        PeptideIndexer -in ${id_file} \\
            -out ${Sample}_${Condition}_${id}_idx.idXML \\
            -threads ${task.cpus} \\
            -fasta ${fasta_decoy} \\
            -decoy_string DECOY \\
            -enzyme:specificity none
    """
}

