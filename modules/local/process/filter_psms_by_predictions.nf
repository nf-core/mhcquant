// Import generic module functions
include { initOptions; saveFiles } from './functions'

params.options = [:]

process FILTER_PSMS_BY_PREDICTIONS {
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:'Intermediate_Results', publish_id:'Intermediate_Results') }

    conda (params.enable_conda ? "bioconda::openms-thirdparty=2.5.0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/openms-thirdparty:2.5.0--6"
    } else {
        container "quay.io/biocontainers/openms-thirdparty:2.5.0--6"
    }
    
    input:
        tuple val(id), val(Sample), val(Condition), file(id_file_psm_filtered)
        tuple val(id), val(Sample), file(peptide_filter_file)

    output:
        tuple val("$id"), val("$Sample"), file("${Sample}_pred_filtered.idXML"), emit: idxml   
        path  "*.version.txt", emit: version

    script:
    """
        IDFilter -in ${id_file_psm_filtered} \\
            -out ${Sample}_pred_filtered.idXML \\
            -whitelist:ignore_modifications \\
            -whitelist:peptides ${peptide_filter_file}\\
            -threads ${task.cpus}

        FileInfo --help &> openms.version.txt
    """
}