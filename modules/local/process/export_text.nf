// Import generic module functions
include { initOptions; saveFiles } from './functions'

params.options = [:]

process EXPORT_TEXT {
    publishDir "${params.outdir}/"

    input:
    tuple val(Sample), file(consensus_resolved) 

    output:
    tuple val(Sample), file("${Sample}.csv") 

    script:
    """
        TextExporter -in ${consensus_resolved} \\
            -out ${Sample}.csv \\
            -threads ${task.cpus} \\
            -id:add_hit_metavalues 0 \\
            -id:add_metavalues 0 \\
            -id:peptides_only
    """
}