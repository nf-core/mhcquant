// Import generic module functions
include { initOptions; saveFiles } from './functions'

params.options = [:]

process EXPORT_MZTAB {
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
        tuple val(Sample), file(feature_file_2)

    output:
        tuple val("id"), val("$Sample"), file("${Sample}.mzTab")

    script:
    """
        MzTabExporter -in ${feature_file_2} \\
            -out ${Sample}.mzTab \\
            -threads ${task.cpus}
    """
}
