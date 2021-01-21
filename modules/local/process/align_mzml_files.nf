// Import generic module functions
include { initOptions; saveFiles } from './functions'

params.options = [:]

process ALIGN_MZML_FILES {

    input:
    tuple val(id), val(Sample), val(Condition), file(mzml_file_align), file(id_file_trafo)

    output:
    tuple val("$id"), val("$Sample"), val("$Condition"), file("${Sample}_${Condition}_${id}_aligned.mzML")

    when:
    !params.skip_quantification

    script:
    """
    MapRTTransformer -in ${mzml_file_align} \\
    -trafo_in ${id_file_trafo} \\
    -out ${Sample}_${Condition}_${id}_aligned.mzML \\
    -threads ${task.cpus}
    """
}