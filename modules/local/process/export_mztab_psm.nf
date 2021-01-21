// Import generic module functions
include { initOptions; saveFiles } from './functions'

params.options = [:]

process EXPORT_MZTAB_PSM {
    publishDir "${params.outdir}/Intermediate_Results/"

    input:
    tuple val(id), val(Sample), val(Condition), file(psm_mztab) 

    output:
    tuple val("$id"), val("$Sample"), file("${Sample}_all_ids_merged.mzTab")

    when:
    params.refine_fdr_on_predicted_subset

    script:
    """
        MzTabExporter -in ${psm_mztab} \\
            -out ${Sample}_all_ids_merged.mzTab \\
            -threads ${task.cpus}
    """

}