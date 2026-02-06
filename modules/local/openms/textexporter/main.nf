process OPENMS_TEXTEXPORTER {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/openms:3.5.0--h78fb946_0' :
        'biocontainers/openms:3.5.0--h78fb946_0' }"

    input:
    tuple val(meta), path(file)

    output:
    tuple val(meta), path("*.tsv"), emit: tsv
    tuple val("${task.process}"), val('openms'), eval("FileInfo --help 2>&1 | sed -nE 's/^Version: ([0-9.]+).*/\\1/p'"), topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args  ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    TextExporter \\
        -in $file \\
        -out ${prefix}_exported.tsv \\
        -threads $task.cpus \\
        -id:add_hit_metavalues 0 \\
        -id:peptides_only \\
        $args
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch ${prefix}.tsv
    """
}
