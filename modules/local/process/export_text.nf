// Import generic module functions
include { initOptions; saveFiles } from './functions'

params.options = [:]

process EXPORT_TEXT {
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