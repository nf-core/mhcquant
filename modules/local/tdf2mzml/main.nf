process TDF2MZML {
    tag "$meta.id"

    container "docker.io/mfreitas/tdf2mzml"

    input:
    tuple val(meta), path(tdf)

    output:
    tuple val(meta), path("*.mzML"), emit: mzml
    tuple val("${task.process}"), val('python'), eval("python3 --version | cut -d ' ' -f2"), topic: versions
    tuple val("${task.process}"), val('tdf2mzml'), eval("echo 0.3.0"), topic: versions

    script:
    def prefix = task.ext.prefix ?: "${tdf.simpleName}"

    """
    tdf2mzml.py -i $tdf -o ${prefix}.mzML
    """

    stub:
    def prefix = task.ext.prefix ?: "${tdf.simpleName}"

    """
    touch ${prefix}.mzML
    """
}
