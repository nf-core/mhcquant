// Import generic module functions
include { initOptions; saveFiles } from './functions'

params.options = [:]

process EXTRACT_PSM_FEATURES_FOR_PERCOLATOR {
    publishDir "${params.outdir}/Intermediate_Results/"

    input:
    tuple val(id), val(Sample), val(Condition), file(id_file_merged)

    output:
    tuple val("$id"), val("$Sample"), val("$Condition"), file("${Sample}_all_ids_merged_psm.idXML")

    script:
    """
    PSMFeatureExtractor -in ${id_file_merged} \\
    -out ${Sample}_all_ids_merged_psm.idXML \\
    -threads ${task.cpus} 
    """

}