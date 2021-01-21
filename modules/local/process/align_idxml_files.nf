// Import generic module functions
include { initOptions; saveFiles } from './functions'

params.options = [:]

process ALIGN_IDXML_FILES {

    input:
    tuple val(id), val(Sample), val(Condition), file(idxml_file_align), file(idxml_file_trafo)

    output:
    tuple val("$id"), val("$Sample"), val("$Condition"), file("${Sample}_${Condition}_${id}_idx_aligned.idXML")

    when:
    !params.skip_quantification

    script:
    """
    MapRTTransformer -in ${idxml_file_align} \\
    -trafo_in ${idxml_file_trafo} \\
    -out ${Sample}_${Condition}_${id}_idx_aligned.idXML \\
    -threads ${task.cpus}
    """
}