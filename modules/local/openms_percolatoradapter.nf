// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process OPENMS_PERCOLATORADAPTER {
    tag "$meta.id"
    label 'process_high'

    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:'Intermediate_Results', publish_id:'Intermediate_Results') }

    conda (params.enable_conda ? "bioconda::openms-thirdparty=2.5.0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/openms-thirdparty:2.5.0--6"
    } else {
        container "quay.io/biocontainers/openms-thirdparty:2.5.0--6"
    }

    input:
        tuple val(meta), path(psm)

    output:
        tuple val(meta), path("*.idXML"), emit: idxml
        path  "*.version.txt", emit: version

    script:
        def software = getSoftwareName(task.process)
        def prefix = options.suffix ? "${meta.id}_${options.suffix}" : "${meta.id}"

        """
            OMP_NUM_THREADS=${task.cpus} \\
            PercolatorAdapter -in ${psm} \\
                -out ${prefix}.idXML \\
                $options.args
            echo \$(FileInfo --help 2>&1) | sed 's/^.*Version: //; s/ .*\$//' &> ${software}-thirdparty.version.txt
        """
}
