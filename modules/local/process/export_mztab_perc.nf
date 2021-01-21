// Import generic module functions
include { initOptions; saveFiles } from './functions'

params.options = [:]

process EXPORT_MZTAB_PERC {
    publishDir "${params.outdir}/Intermediate_Results/"

    input:
    tuple val(id), val(Sample), file(percolator_mztab)

    output:
    tuple val("$id"), val("$Sample"), file("${Sample}_all_ids_merged_psm_perc_filtered.mzTab")

    when:
    params.refine_fdr_on_predicted_subset

    script:
    """
        MzTabExporter -in ${percolator_mztab} \\
            -out ${Sample}_all_ids_merged_psm_perc_filtered.mzTab \\
            -threads ${task.cpus}
    """
// 
}