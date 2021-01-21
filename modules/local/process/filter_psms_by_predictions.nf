// Import generic module functions
include { initOptions; saveFiles } from './functions'

params.options = [:]

process FILTER_PSMS_BY_PREDICTIONS {
    publishDir "${params.outdir}/Intermediate_Results/"
    
    input:
    tuple val(id), val(Sample), val(Condition), file(id_file_psm_filtered)
    tuple val(id), val(Sample), file(peptide_filter_file)

    output:
    tuple val("$id"), val("$Sample"), file("${Sample}_pred_filtered.idXML")

    when:
    params.refine_fdr_on_predicted_subset    

    script:
    """
        IDFilter -in ${id_file_psm_filtered} \\
            -out ${Sample}_pred_filtered.idXML \\
            -whitelist:ignore_modifications \\
            -whitelist:peptides ${peptide_filter_file}\\
            -threads ${task.cpus} \\
    """
}