// Import generic module functions
include { initOptions; saveFiles } from './functions'

params.options = [:]

process ALIGN_IDS {

    input:
    tuple val(id), val(Sample), val(Condition), file(id_names)

    output:
    tuple val("$Sample"), file("*.trafoXML")

    when:
    !params.skip_quantification

    script:
    def out_names = id_names.collect { it.baseName+'.trafoXML' }.join(' ')

    """
    MapAlignerIdentification -in $id_names \\
    -trafo_out $out_names \\
    -model:type linear \\
    -algorithm:max_rt_shift ${params.max_rt_alignment_shift}
    """
}