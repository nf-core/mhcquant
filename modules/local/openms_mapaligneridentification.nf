// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options    = initOptions(params.options)

process OPENMS_MAPALIGNERIDENTIFICATION {

    conda (params.enable_conda ? "bioconda::openms=2.5.0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/openms:2.5.0--h4afb90d_6"
    } else {
        container "quay.io/biocontainers/openms:2.5.0--h4afb90d_6"
    }

    input:
        tuple val(id), val(Sample), val(Condition), path(id_names)

    output:
        tuple val("$Sample"), path("*.trafoXML"), emit: trafoxml   
        path  "*.version.txt", emit: version

    script:
        def out_names = id_names.collect { it.baseName+'.trafoXML' }.join(' ')
        def software = getSoftwareName(task.process)

        """
            MapAlignerIdentification -in $id_names \\
                -trafo_out $out_names \\
                $options.args 
            echo \$(FileInfo --help 2>&1) | sed 's/^.*Version: //; s/ .*\$//' &> ${software}.version.txt
        """
}