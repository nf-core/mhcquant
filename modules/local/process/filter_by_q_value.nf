// Import generic module functions
include { initOptions; saveFiles } from './functions'

params.options = [:]

process FILTER_BY_Q_VALUE {
    publishDir "${params.outdir}/Intermediate_Results/"
    
    input:
    tuple val(id), val(Sample), val(Condition), file(id_file_perc) 

    output:
    tuple val("$id"), val("$Sample"), file("${Sample}_all_ids_merged_psm_perc_filtered.idXML")

    // when:
    // !params.refine_fdr_on_predicted_subset

    script:
    """
        IDFilter -in ${id_file_perc} \\
            -out ${Sample}_all_ids_merged_psm_perc_filtered.idXML \\
            -threads ${task.cpus} \\
            -score:pep ${params.fdr_threshold} \\
            -remove_decoys \\
            -precursor:length '${params.peptide_min_length}:${params.peptide_max_length}' \\
            -delete_unreferenced_peptide_hits
    """
}