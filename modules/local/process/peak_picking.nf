// Import generic module functions
include {  saveFiles } from './functions'

params.options = [:]

process PEAK_PICKING {

    input:
    tuple val(id), val(Sample), val(Condition), file(mzml_unpicked) //from input_mzmls_unpicked

    output:
    tuple val("$id"), val("$Sample"), val("$Condition"), file("${mzml_unpicked.baseName}.mzML") //into (input_mzmls_picked, input_mzmls_align_picked)

    when:
    params.run_centroidisation

    script:
     """
     PeakPickerHiRes -in ${mzml_unpicked} \\
                     -out ${mzml_unpicked.baseName}.mzML \\
                     -algorithm:ms_levels ${params.pick_ms_levels}
     """
}