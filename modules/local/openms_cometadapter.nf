// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process OPENMS_COMETADAPTER {
    tag "$meta.id"
    label 'process_high'

    conda (params.enable_conda ? "bioconda::openms-thirdparty=2.5.0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/openms-thirdparty:2.5.0--6"
    } else {
        container "quay.io/biocontainers/openms-thirdparty:2.5.0--6"
    }

    input:
        tuple val(meta), path(mzml), path(fasta)

    output:
        tuple val(meta), path("*.idXML"), emit: idxml
        path  "*.version.txt", emit: version

    script:
        def software = getSoftwareName(task.process)
        def prefix   = options.suffix ? "${mzml.baseName}_${options.suffix}" : "${mzml.baseName}"

        """
            CometAdapter -in ${mzml} \\
                -out ${prefix}.idXML \\
                -database ${fasta} \\
                -threads $task.cpus $options.args
            echo \$(FileInfo --help 2>&1) | sed 's/^.*Version: //; s/ .*\$//' &> ${software}-thirdparty.version.txt
        """
}
