// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process OPENMS_MAPALIGNERIDENTIFICATION {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::openms=2.5.0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/openms:2.5.0--h4afb90d_6"
    } else {
        container "quay.io/biocontainers/openms:2.5.0--h4afb90d_6"
    }

    input:
        tuple val(meta), path(idxml)

    output:
        tuple val(meta), path("*.trafoXML"), emit: trafoxml   
        path  "*.version.txt", emit: version

    script:
        def software = getSoftwareName(task.process)
        def out_names = idxml.collect { it.baseName+'.trafoXML' }.join(' ')

        """
            MapAlignerIdentification -in ${idxml} \\
                -trafo_out ${out_names} \\
                $options.args 
            echo \$(FileInfo --help 2>&1) | sed 's/^.*Version: //; s/ .*\$//' &> ${software}.version.txt
        """
}