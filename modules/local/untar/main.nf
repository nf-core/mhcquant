process UNTAR {
    tag "$archive"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
        'nf-core/ubuntu:20.04' }"

    input:
    tuple val(meta), path(archive)

    output:
    tuple val(meta), path("*.d"), emit: untar
    path "versions.yml"         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args  = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: archive.baseName.replaceAll(/\.tar(\.gz)?$/, '')

    """
    mkdir $prefix
    depth=\$(tar -tf "${archive}" | grep '\\.d/\$' | head -n 1 | tr -cd '/' | wc -c)

    tar \\
        -C $prefix \\
        -xavf \\
        $args \\
        $archive \\
        --strip-components=\$depth \\
        $args2

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        untar: \$(echo \$(tar --version 2>&1) | sed 's/^.*(GNU tar) //; s/ Copyright.*\$//')
    END_VERSIONS
    """

    stub:
    prefix    = task.ext.prefix ?: ( meta.id ? "${meta.id}" : archive.toString().replaceFirst(/\.[^\.]+(.gz)?$/, ""))
    """
    mkdir $prefix
    touch ${prefix}/file.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        untar: \$(echo \$(tar --version 2>&1) | sed 's/^.*(GNU tar) //; s/ Copyright.*\$//')
    END_VERSIONS
    """
}
