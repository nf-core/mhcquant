// Import generic module functions
include { initOptions; saveFiles } from './functions'

params.options = [:]

//TODO: combine in a subflow --> when needs to be removed
process EXPORT_MZTAB_PERC {
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
        tuple val(id), val(Sample), file(percolator_mztab)

    output:
        tuple val("$id"), val("$Sample"), file("${Sample}_all_ids_merged_psm_perc_filtered.mzTab"), emit: mztab   
        path  "*.version.txt", emit: version

    when:
        params.refine_fdr_on_predicted_subset

    script:
    """
        MzTabExporter -in ${percolator_mztab} \\
            -out ${Sample}_all_ids_merged_psm_perc_filtered.mzTab \\
            -threads ${task.cpus}

        FileInfo --help &> openms.version.txt
    """
}