// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]

process OPENMS_IDCONFLICTRESOLVER {
    conda (params.enable_conda ? "bioconda::openms=2.5.0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/openms:2.5.0--h4afb90d_6"
    } else {
        container "quay.io/biocontainers/openms:2.5.0--h4afb90d_6"
    }
    
    input:
        tuple val(Sample), path(consensus)

    output:
        tuple val(Sample), path("${Sample}_resolved.consensusXML"), emit: consensusxml   
        path  "*.version.txt", emit: version
        
    script:
        def software = getSoftwareName(task.process)
        
        """
            IDConflictResolver -in ${consensus} \\
                -out ${Sample}_resolved.consensusXML \\
                -threads ${task.cpus}
            echo \$(FileInfo --help 2>&1) | sed 's/^.*Version: //; s/ .*\$//' &> ${software}.version.txt
        """
}