// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process OPENMS_MAPRTTRANSFORMER {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::openms=2.5.0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/openms:2.5.0--h4afb90d_6"
    } else {
        container "quay.io/biocontainers/openms:2.5.0--h4afb90d_6"
    }

    input:
        tuple val(meta), path(alignment_file), path(trafoxml)

    output:   
        tuple val(meta), path("*_aligned.*"), emit: aligned   
        path  "*.version.txt", emit: version

    script:
        def software = getSoftwareName(task.process)
        def fileExt = alignment_file.collect { it.name.tokenize("\\.")[1] }.join(' ')

        """
            MapRTTransformer -in ${alignment_file} \\
                -trafo_in ${trafoxml} \\
                -out ${meta.id}_aligned.${fileExt} \\
                -threads $task.cpus
            echo \$(FileInfo --help 2>&1) | sed 's/^.*Version: //; s/ .*\$//' &> ${software}.version.txt
        """
}