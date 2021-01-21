// Import generic module functions
include { initOptions; saveFiles } from './functions'

params.options = [:]

process RESOLVE_CONFLICTS {
    
    input:
    tuple val(Sample), file(consensus)

    output:
    tuple val(Sample), file("${Sample}_resolved.consensusXML")

    when:
    !params.skip_quantification

    script:
    """
        IDConflictResolver -in ${consensus} \\
            -out ${Sample}_resolved.consensusXML \\
            -threads ${task.cpus}
     """

}