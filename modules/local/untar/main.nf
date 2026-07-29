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
    tuple val("${task.process}"), val('untar'), eval("tar --version 2>&1 | grep -oE -m1 '[0-9]+\\.[0-9.]+'"), topic: versions

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
    """

    stub:
    prefix    = task.ext.prefix ?: ( meta.id ? "${meta.id}" : archive.toString().replaceFirst(/\.[^\.]+(.gz)?$/, ""))

    """
    mkdir $prefix
    touch ${prefix}/file.txt
    """
}
