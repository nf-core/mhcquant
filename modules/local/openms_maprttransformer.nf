process OPENMS_MAPRTTRANSFORMER {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::openms=2.6.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/openms:2.6.0--h4afb90d_0' :
        'quay.io/biocontainers/openms:2.6.0--h4afb90d_0' }"

    input:
        tuple val(meta), path(alignment_file), path(trafoxml)

    output:
        tuple val(meta), path("*_aligned.*"), emit: aligned
        path "versions.yml"                 , emit: versions

    script:
        def fileExt          = alignment_file.collect { it.name.tokenize("\\.")[1] }.join(' ')

        """
        MapRTTransformer -in $alignment_file \\
            -trafo_in $trafoxml \\
            -out ${meta.id}_aligned.${fileExt} \\
            -threads $task.cpus

        cat <<-END_VERSIONS > versions.yml
        ${task.process}:
            openms: \$(echo \$(FileInfo --help 2>&1) | sed 's/^.*Version: //; s/-.*\$//' | sed 's/ -*//; s/ .*\$//')
        END_VERSIONS
        """
}
