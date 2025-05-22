process OPENMS_MAPRTTRANSFORMER {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/openms:3.4.0--hc77a4c7_0' :
        'biocontainers/openms:3.4.0--hc77a4c7_0' }"

    input:
    tuple val(meta), path(alignment_file), path(trafoxml)

    output:
    tuple val(meta), path("*_aligned.*"), emit: aligned
    path "versions.yml"                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args    = task.ext.args  ?: ''
    def prefix  = task.ext.prefix ?: "${meta.id}"
    def fileExt = alignment_file.collect { it.name.tokenize("\\.")[1] }.join(' ')

    """
    MapRTTransformer \\
        -in $alignment_file \\
        -trafo_in $trafoxml \\
        -out ${prefix}.${fileExt} \\
        -threads $task.cpus \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        openms: \$(echo \$(FileInfo --help 2>&1) | sed 's/^.*Version: //; s/-.*\$//' | sed 's/ -*//; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    def prefix  = task.ext.prefix ?: "${meta.id}"
    def fileExt = alignment_file.collect { it.name.tokenize("\\.")[1] }.join(' ')

    """
    touch ${prefix}.${fileExt}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        openms: \$(echo \$(FileInfo --help 2>&1) | sed 's/^.*Version: //; s/-.*\$//' | sed 's/ -*//; s/ .*\$//')
    END_VERSIONS
    """
}
