// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options    = initOptions(params.options)

process OPENMS_PERCOLATORADAPTER {
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
        tuple val(id), val(Sample), val(Condition), path(psm_file)

    output:
        tuple val("$id"), val("$Sample"), val("$Condition"), path("*.idXML"), emit: idxml   
        path  "*.version.txt", emit: version

    script:
        def software = getSoftwareName(task.process)
        def prefix = options.suffix ? "${Sample}_${options.suffix}" : "${Sample}_${id}"

        """
            OMP_NUM_THREADS=${task.cpus} \\
            PercolatorAdapter -in ${psm_file} \\
                -out ${prefix}.idXML \\
                $options.args
            echo \$(FileInfo --help 2>&1) | sed 's/^.*Version: //; s/ .*\$//' &> ${software}.version.txt
        """
}