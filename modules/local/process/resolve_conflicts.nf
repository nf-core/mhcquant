// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]

//TODO: combine in a subflow --> when needs to be removed
process RESOLVE_CONFLICTS {
    conda (params.enable_conda ? "bioconda::openms-thirdparty=2.5.0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/openms-thirdparty:2.5.0--6"
    } else {
        container "quay.io/biocontainers/openms-thirdparty:2.5.0--6"
    }
    
    input:
        tuple val(Sample), file(consensus)

    output:
        tuple val(Sample), file("${Sample}_resolved.consensusXML"), emit: consensusxml   
        path  "*.version.txt", emit: version

    when:
        !params.skip_quantification

    script:
        """
        IDConflictResolver -in ${consensus} \\
                -out ${Sample}_resolved.consensusXML \\
                -threads ${task.cpus}

        FileInfo --help &> openms.version.txt
        """
}