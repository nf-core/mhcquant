// Import generic module functions
include { initOptions; saveFiles } from './functions'

params.options = [:]

process CALCULATE_FDR_FOR_ID_ALIGNMENT  {
    input:
    tuple val(id), val(Sample), val(Condition), file(id_file_idx)

    output:
    tuple val("$id"), val("$Sample"), val("$Condition"), file("${Sample}_${Condition}_${id}_idx_fdr.idXML") 

    when:
    !params.skip_quantification

    script:
    """
    FalseDiscoveryRate -in ${id_file_idx} \\
    -protein 'false' \\
    -out ${Sample}_${Condition}_${id}_idx_fdr.idXML \\
    -threads ${task.cpus}
    """
}