// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process OPENMS_RTMODEL {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::openms=2.5.0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/openms:2.5.0--h4afb90d_6"
    } else {
        container "quay.io/biocontainers/openms:2.5.0--h4afb90d_6"
    }

    input:
        tuple val(meta), path(rt_training)

    output:
        tuple val(meta), path("*_rt_training.txt"), path("*.paramXML"), path("*_trainset.txt"), emit: complete
        path  "*.version.txt", emit: version

    script:
        def software = getSoftwareName(task.process)
        def prefix = options.suffix ? "${meta.sample}_${options.suffix}" : "${meta.sample}"

        """
            RTModel -in ${rt_training} \\
                -cv:skip_cv \\
                -out ${prefix}_rt_training.txt \\
                -out_oligo_params ${prefix}_params.paramXML \\
                -out_oligo_trainset ${prefix}_trainset.txt
            echo \$(FileInfo --help 2>&1) | sed 's/^.*Version: //; s/ .*\$//' &> ${software}.version.txt
        """
}
