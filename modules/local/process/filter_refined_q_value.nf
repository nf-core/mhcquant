// Import generic module functions
include { initOptions; saveFiles } from './functions'

params.options = [:]

//TODO: combine in a subflow --> when needs to be removed
process FILTER_REFINED_Q_VALUE {
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
        tuple val(id), val(Sample), file(id_file_perc_pred) 
    
    output:
        tuple val("$id"), val("$Sample"), file("${Sample}_perc_subset_filtered.idXML"), emit: idxml   
        path  "*.version.txt", emit: version

    when:
        params.refine_fdr_on_predicted_subset     

    script:
    """      
        IDFilter -in ${id_file_perc_pred} \\
            -out ${Sample}_perc_subset_filtered.idXML \\
            -threads ${task.cpus} \\
            -score:pep ${params.fdr_threshold} \\
            -remove_decoys \\
            -precursor:length '${params.peptide_min_length}:${params.peptide_max_length}' \\
            -delete_unreferenced_peptide_hits

        FileInfo --help &> openms.version.txt
    """

}