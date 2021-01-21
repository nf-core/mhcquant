
// Import generic module functions
include { initOptions; saveFiles } from './functions'

params.options = [:]

process EXPORT_MZTAB {
    publishDir "${params.outdir}/"

    input:
    tuple val(Sample), file(feature_file_2)

    output:
    tuple val("id"), val("$Sample"), file("${Sample}.mzTab")

    script:
    """
        MzTabExporter -in ${feature_file_2} \\
            -out ${Sample}.mzTab \\
            -threads ${task.cpus}
    """
}
