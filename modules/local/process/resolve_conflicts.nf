// Import generic module functions
include { initOptions; saveFiles } from './functions'

params.options = [:]

process RESOLVE_CONFLICTS {
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:'') }

    conda (params.enable_conda ? "bioconda::openms-thirdparty=2.5.0" : null)
    
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/openms-thirdparty:2.5.0--6"
    } else {
        container "quay.io/biocontainers/openms-thirdparty:2.5.0--6"
    }
    
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